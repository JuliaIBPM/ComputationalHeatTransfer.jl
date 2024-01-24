module ComputationalHeatTransfer

# using DocStringExtensions
using Reexport
using UnPack
using ImmersedLayers

# reexport exports them further from these packages
@reexport using ImmersedLayers
 
export DirichletHeatConductionProblem
export timestep_fourier
export heatconduction_ode_rhs
export heatconduction_bc_rhs
export heatconduction_constraint_force
export heatconduction_bc_op

export heatconduction_rhs
export timestep_fourier_cfl
export UnboundedHeatConductionProblem

export NeumannHeatConductionProblem
export heatconduction_ode_implicit_rhs
export heatconduction_ode_explicit_rhs
export heatconduction_neumann_bc_rhs
export heatconduction_neumann_constraint_force
export heatconduction_neumann_bc_op
export heatconduction_neumann_bc_reg

function default_convection_velocity(v,t,base_cache,phys_params)
    fill!(v,0)
end

const DEFAULT_CONVECTION_VELOCITY_FUNCTION = default_convection_velocity

function get_convection_velocity_function(phys_params::Dict)
    return get(phys_params,"convection velocity model",DEFAULT_CONVECTION_VELOCITY_FUNCTION)
end

function get_heating_model(forcing::Dict)
	return get(forcing,"heating models",nothing)
end

function get_heating_model(forcing::Nothing)
	return get_heating_model(Dict())
end

#####################################################################

# Calculate the RHS of the ODE
function heatconduction_ode_rhs!(dT,T,x,sys::ILMSystem,t)
    @unpack bc, forcing, phys_params, extra_cache, base_cache = sys
    @unpack dTb, dT_tmp, v, cdcache, fcache = extra_cache

    κ = phys_params["diffusivity"]

    # Calculate the double-layer term
    prescribed_surface_jump!(dTb,x,t,sys)
    surface_divergence!(dT,-κ*dTb,sys)

    # This provides the convection velocity at time `t`
    convection_function = get_convection_velocity_function(phys_params)
    convection_function(v,t,base_cache,phys_params)

    # This provides the heating model at time `t`
    # heating_function = get_heating_model_function(phys_params)
    # heating_function(T,x,fcache,phys_params)

    # Compute the convective derivative term `N(v,T)`
    fill!(dT_tmp,0.0)
    convective_derivative!(dT_tmp,v,T,base_cache,cdcache)
    dT .-= dT_tmp

    # Compute the contribution from the forcing models to the right-hand side
    fill!(dT_tmp,0.0)
    apply_forcing!(dT_tmp,T,x,t,fcache,sys)
    dT .+= dT_tmp

    return dT
end

# Calculate the RHS of the boundary condition
function heatconduction_bc_rhs!(dTb,x,sys::ILMSystem,t)
    prescribed_surface_average!(dTb,x,t,sys)
    return dTb
end

# Calculate the contribution to dT/dt
function heatconduction_constraint_force!(dT,σ,x,sys::ILMSystem)
    regularize!(dT,-σ,sys)
    return dT
end

# Interpolate temperature vector onto the boundary
function heatconduction_bc_op!(dTb,T,x,sys::ILMSystem)
    interpolate!(dTb,T,sys)
    return dTb
end

# Timestep function for when convective time is irrelevant
function timestep_fourier(u,sys)
    @unpack phys_params = sys
    g = get_grid(sys)
    κ = phys_params["diffusivity"]
    Fo = phys_params["Fourier"]
    Δt = Fo*cellsize(g)^2/κ
    return Δt
end

# Timestep function for when convective time is relevant
function timestep_fourier_cfl(u,sys)
    @unpack phys_params = sys
    g = get_grid(sys)
    κ = phys_params["diffusivity"]
    Ω = phys_params["angular velocity"]
    L = phys_params["length scale"]
    Fo = phys_params["Fourier"]
    Co = phys_params["CFL"]
    Δt = min(Fo*cellsize(g)^2/κ,Co*cellsize(g)/Ω/L)
    return Δt
end

# Set up extra cache - Dirichlet problem
@ilmproblem DirichletHeatConduction scalar

# Structure cache for DirichletHeatConduction
struct DirichletHeatConductionCache{VT,CDT,FRT,DTST,DTT,FT} <: AbstractExtraILMCache
    v :: VT
    cdcache :: CDT
    fcache :: FRT
    dTb :: DTST
    dT_tmp :: DTT
    f :: FT
end

# Testing Method Dispatch - Dirichlet and Unbounded
function ImmersedLayers.prob_cache(prob::DirichletHeatConductionProblem,
                                    base_cache::BasicILMCache{N,scaling}) where {N,scaling}
    @unpack phys_params, forcing = prob
    @unpack gdata_cache, g = base_cache

    dTb = zeros_surface(base_cache)

    # Construct a Laplacian outfitted with the diffusivity
    κ = phys_params["diffusivity"]
    heat_L = Laplacian(base_cache,κ)

    # Create cache for the convective derivative
    v = zeros_gridgrad(base_cache)
    cdcache = ConvectiveDerivativeCache(base_cache)

    # Create cache for the forcing regions
    heating_model = get_heating_model(forcing)
    fcache = ForcingModelAndRegion(heating_model,base_cache)

    # The state here is temperature
    f = _get_dirichlet_ode_function_list(heat_L,base_cache)

    dT_tmp = zeros_grid(base_cache)

    DirichletHeatConductionCache(v,cdcache,fcache,dTb,dT_tmp,f)

end

_get_dirichlet_ode_function_list(heat_L,base_cache::BasicILMCache{N}) where {N} =
            ODEFunctionList(state = zeros_grid(base_cache),
                            constraint = zeros_surface(base_cache),
                            ode_rhs = heatconduction_ode_rhs!,
                            lin_op = heat_L,
                            bc_rhs = heatconduction_bc_rhs!,
                            constraint_force = heatconduction_constraint_force!,
                            bc_op = heatconduction_bc_op!)

_get_dirichlet_ode_function_list(heat_L,base_cache::BasicILMCache{0}) =
            ODEFunctionList(state = zeros_grid(base_cache),
                            ode_rhs = heatconduction_ode_rhs!,
                            lin_op = heat_L)

#####################################################################

# Implementation of R[q] the implicit part of the equation RHS
function heatconduction_ode_implicit_rhs!(dT,x,sys::ILMSystem,t)
    @unpack bc, forcing, phys_params, extra_cache, base_cache = sys
    @unpack dqbtmp = extra_cache

    fill!(dT,0.0)
    fill!(dqbtmp,0.0)

    # Calculate the single-layer term on the RHS
    prescribed_surface_jump!(dqbtmp,x,t,sys)
    regularize!(dT,dqbtmp,sys)

    return dT
end

# Implementation of explicit part of equation RHS
function heatconduction_ode_explicit_rhs!(dT,T,x,sys::ILMSystem,t)
    @unpack extra_cache, base_cache, phys_params, forcing = sys
    @unpack v, cdcache, dT_tmp, fcache, Ttmp = extra_cache

    fill!(dT,0.0)
    # Compute the contribution from the forcing models to the right-hand side
    apply_forcing!(dT,T,x,t,fcache,sys)

    # This provides the convection velocity at time `t`
    convection_function = get_convection_velocity_function(phys_params)
    convection_function(v,t,base_cache,phys_params)

    # Compute the convective derivative term `N(v,T)`
    fill!(dT_tmp,0.0)
    convective_derivative!(dT_tmp,v,T,base_cache,cdcache)
    dT .-= dT_tmp

    return dT
end

# Calculate RHS of the boundary condition
function heatconduction_neumann_bc_rhs!(dqb,x,sys::ILMSystem,t)
    fill!(dqb,0.0)
    prescribed_surface_average!(dqb,x,t,sys)
    return dqb
end

# Calculates negative contribution to dT/dt from Lagrange multiplier
function heatconduction_neumann_constraint_force!(dT,Tjump,x,sys::ILMSystem)
    @unpack phys_params = sys
    κ = phys_params["diffusivity"]
    fill!(dT,0.0)
    surface_divergence!(dT,κ*Tjump,sys)
    return dT
end

# Calculate LHS transpose - surface gradient operation
function heatconduction_neumann_bc_op!(dqb,T,x,sys::ILMSystem)
    @unpack phys_params = sys
    κ = phys_params["diffusivity"]
    fill!(dqb,0.0)
    surface_grad!(dqb,T,sys)
    dqb .*= -κ
    return dqb
end

# Compute the final left-hand-side operation of constraint equations
function heatconduction_neumann_bc_reg!(dqb,Tjump,x,sys::ILMSystem)
    @unpack extra_cache, phys_params = sys
    @unpack Ttmp, vtmp = extra_cache
    κ = phys_params["diffusivity"]

    fill!(vtmp,0.0)
    fill!(dqb,0.0)
    regularize_normal!(vtmp,κ*Tjump,sys)
    normal_interpolate!(dqb,vtmp,sys)

    return dqb
end

# Set up extra cache - Neumann problem
@ilmproblem NeumannHeatConduction scalar

# Structure cache for NeumannHeatConduction
struct NeumannHeatConductionCache{VT,CDT,DTT,GT,GVT,FRT,FT} <: AbstractExtraILMCache
    v :: VT
    cdcache :: CDT
    dT_tmp :: GT
    Ttmp :: GT
    vtmp :: GVT
    dqbtmp :: DTT
    fcache :: FRT
    f :: FT
end

# Method dispatch - problem cache for NeumannHeatConductionProblem
function ImmersedLayers.prob_cache(prob::NeumannHeatConductionProblem,
                                base_cache::BasicILMCache{N,scaling}) where {N,scaling}

    @unpack phys_params, forcing = prob
    @unpack gdata_cache, g = base_cache

    Ttmp = zeros_grid(base_cache)
    vtmp = zeros_gridgrad(base_cache)
    dqbtmp = zeros_surface(base_cache)

    # Construct a Laplacian outfitted with the diffusivity
    κ = phys_params["diffusivity"]
    heat_L = Laplacian(base_cache,κ)

    # Create cache for the convective derivative
    v = zeros_gridgrad(base_cache)
    cdcache = ConvectiveDerivativeCache(base_cache)

    # Create cache for the forcing regions
    fcache = ForcingModelAndRegion(forcing["heating models"],base_cache)

    dT_tmp = zeros_grid(base_cache)

    # State (grid temperature data) and constraint (surface Lagrange multipliers)
    f = _get_neumann_ode_function_list(heat_L,base_cache)

    NeumannHeatConductionCache(v,cdcache,dT_tmp,Ttmp,vtmp,dqbtmp,fcache,f)
end

_get_neumann_ode_function_list(heat_L,base_cache::BasicILMCache{N}) where {N} =
                ODEFunctionList(state = zeros_grid(base_cache),
                                constraint = zeros_surface(base_cache),
                                ode_rhs = heatconduction_ode_explicit_rhs!,
                                lin_op = heat_L,
                                bc_rhs = heatconduction_neumann_bc_rhs!,
                                constraint_force = heatconduction_neumann_constraint_force!,
                                bc_op = heatconduction_neumann_bc_op!,
                                bc_regulator = heatconduction_neumann_bc_reg!,
                                ode_implicit_rhs = heatconduction_ode_implicit_rhs!)
     
_get_neumann_ode_function_list(heat_L,base_cache::BasicILMCache{0}) =
                ODEFunctionList(state = zeros_grid(base_cache),
                                ode_rhs = heatconduction_rhs!,
                                lin_op = heat_L)

end
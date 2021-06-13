#### Operators for a system with a body and stationary points ####

export temperature_bc_constraint_rhs!, heatflux_bc_constraint_rhs!

function temperature_bc_constraint_rhs!(Ts::ScalarData{N},sys::HeatConduction{NX,NY,N},t::Real) where {NX,NY,N}
    surface_temperature!(sys.ΔTs,sys,t)
    Ts .= sys.ΔTs
    return Ts
end


function heatflux_bc_constraint_rhs!(qs::ScalarData{N},sys::HeatConduction{NX,NY,N},t::Real) where {NX,NY,N}
    fill!(qs,0.0)
    # for now, it's adiabatic only
    return qs
end

# Constraint operators, using stored regularization and interpolation operators
# B₁ᵀ = Eᵀ
function temperature_constraint_force!(out::Nodes{Primal,NX,NY},σ::ScalarData{N},sys::HeatConduction{NX,NY,N}) where {NX,NY,N}
    out .= sys.Rc*σ
    return out
end

function temperature_constraint_force!(out::Nodes{Primal,NX,NY},σ,sys::HeatConduction{NX,NY,0}) where {NX,NY}
    return out
end

function heatflux_constraint_force!(out::Nodes{Primal,NX,NY},τ::ScalarData{N},sys::HeatConduction{NX,NY,N}) where {NX,NY,N}
    @unpack params, normals, Rf, Vb, Vf = sys
    @unpack d, α = params
    fact = α/cellsize(sys)
    #fact = isnothing(d) ? fact : fact/d

    fill!(out,0.0)
    product!(Vb,normals,τ)
    Vf .= Rf*Vb
    divergence!(out,Vf)
    out .*= fact
    return out
end

function heatflux_constraint_force!(out::Nodes{Primal,NX,NY},σ,sys::HeatConduction{NX,NY,0}) where {NX,NY}
    return out
end

#=
function _vel_ns_op_constraint_force!(u::Edges{Primal,NX,NY},τ,sys::HeatConduction{NX,NY,0}) where {NX,NY}
    return u
end

function _vel_ns_op_constraint_force!(u::Edges{Primal,NX,NY},τ::VectorData{N},sys::HeatConduction{NX,NY,N}) where {NX,NY,N}
    u .= sys.Rf*τ
    return u
end
=#

# B₂ = E
function temperature_constraint_op!(out::ScalarData{N},T::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY,N}) where {NX,NY,N}
    out .= sys.Ec*T
    return out
end

function heatflux_constraint_op!(out::ScalarData{N},T::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY,N}) where {NX,NY,N}
    @unpack params, normals, Ef, gradTb, gradTf = sys
    @unpack k = params
    fill!(gradTf,0.0)

    grad!(gradTf,T)
    gradTb .= Ef*gradTf
    pointwise_dot!(out,gradTb,normals)
    out .*= -k/cellsize(sys)
    return out
end

import Base: size

using UnPack

import CartesianGrids: cellsize, origin
import RigidBodyTools: assign_velocity!
import ImmersedLayers: normals, areas

"""
    HeatConductionParameters(ρ::Float64,c::Float64,k::Float64,[thickness=nothing])

Set the density `ρ`, specific heat `c`, heat conductivity `k`, and optionally,
the thickness.
"""
struct HeatConductionParameters
    ρ :: Float64
    c :: Float64
    k :: Float64
    α :: Float64
    d :: Union{Float64,Nothing}
end
HeatConductionParameters(ρ,c,k,d) = HeatConductionParameters(ρ,c,k,k/(ρ*c),d)
HeatConductionParameters(ρ,c,k;thickness=nothing) = HeatConductionParameters(ρ,c,k,thickness)


include("operators/heaters.jl")
include("operators/linesource.jl")


"""
$(TYPEDEF)

A system type that utilizes a grid of `NX` x `NY` dual cells and `N` Lagrange forcing
points to solve the discrete transient heat conduction equation. The
parameter `static_points` specifies whether the forcing points remain static in the
grid. It should be set to `false` if a supplied motion requires that the points move.

# Constructors:

`HeatConduction(params,Δx,xlimits,ylimits,Δt
              [,pulses=nothing])` specifies the Reynolds number `Re`, the grid
              spacing `Δx`, the dimensions of the domain in the tuples `xlimits`
              and `ylimits` (excluding the ghost cells), and the time step size `Δt`.
              The other arguments are optional. The `freestream` argument can be
              passed as either a tuple (a static freestream) or a `RigidBodyMotion`
              for a time-varying freestream. The `pulses` argument can be
              used to pass in one or more spatiotemporal pulses.


`HeatConduction(Re,Δx,xlimits,ylimits,Δt,bodies::Body/BodyList
              [,problem_side=ExternalInternalProblem]
              [,ddftype=CartesianGrids.Yang3])` passes the body
              information. This constructor
              sets the motions of the body/ies to be stationary.
              The same optional arguments used for the basic constructor
              also apply for this one. In addition, the `problem_side` can be set to
              `ExternalProblem` (default), `InternalProblem`, or `ExternalInternalProblem`.
              However, it is forced to `ExternalInternalProblem` for open Bodies
              (like `Plate` type).

`HeatConduction(Re,Δx,xlimits,ylimits,Δt,bodies::Body/BodyList,
              motions::RigidBodyMotion/RigidMotionList
              [,static_points=false])`
              passes the body and associated motion information.
              The list of motions must be the same length as the list of bodies.
              The same optional arguments used for the other constructors
              also apply for this one. In addition, `static_points` can
              be set to `true` if the supplied motion should not cause the
              points to move.

"""
mutable struct HeatConduction{NX, NY, N, MT<:PointMotionType, SD<:ProblemSide, DDF<:CartesianGrids.DDFType, FT, SP} #, RKT}
    # Physical Parameters
    "Heat transfer parameters"
    params::HeatConductionParameters

    "Bodies"
    bodies::Union{BodyList,Nothing}
    "BC type"
    bctype::Type
    "Body temperatures"
    bodytemps::Union{RigidMotionList,Nothing}
    "Volumetric heat inputs"
    qforce::Union{Vector{ModulatedField},Nothing}

    "Prescribed area heat flux"
    qflux::Union{Vector{PrescribedHeatFlux},Nothing}
    "Prescribed area h*dT heat model"
    qhdT::Union{Vector{PrescribedHeatModel},Nothing}
    "Prescribed line source"
    qline::Union{Vector{PrescribedLineSource},Nothing}

    # Discretization
    "Grid metadata"
    grid::CartesianGrids.PhysicalGrid{2}
    "Time step"
    Δt::Float64

    # Layers
    dlc::Union{DoubleLayer,Nothing} # used for heat flux surface terms

    # Coordinate data, if present
    points::VectorData{N,Float64}
    normals::VectorData{N,Float64}

    # Pre-stored regularization and interpolation matrices (if present)
    Rc::Union{RegularizationMatrix,Nothing} # cell centers
    Ec::Union{InterpolationMatrix,Nothing}
    Rf::Union{RegularizationMatrix,Nothing} # cell faces
    Ef::Union{InterpolationMatrix,Nothing}
    Cc::Union{Matrix,Nothing}

    # Operators
    f :: FT

    # state vector
    state_prototype :: SP

    # Cache space
    Sc::Nodes{Primal, NX, NY,Float64}
    Sb::ScalarData{N,Float64}
    Vf::Edges{Primal, NX, NY,Float64}
    gradTf::Edges{Primal, NX, NY,Float64}
    Vb::VectorData{N,Float64}
    gradTb::VectorData{N,Float64}
    ΔTs::ScalarData{N,Float64}
    σ::ScalarData{N,Float64}
    τ::ScalarData{N,Float64}

end

function HeatConduction(params::HeatConductionParameters, Δx::Real, xlimits::Tuple{Real,Real},ylimits::Tuple{Real,Real}, Δt::Real;
                       bodies::Union{BodyList,Nothing} = nothing,
                       bodytemps::Union{RigidMotionList,Nothing} = nothing,
                       qforce::PT = nothing,
                       qflux = nothing,
                       qmodel = nothing,
                       qline = nothing,
                       static_points = true,
                       bctype = AdiabaticBC,
                       problem_side::Type{SD} = InternalProblem,
                       ddftype=CartesianGrids.Yang3) where {PT,SD<:ProblemSide}

    @unpack α = params

    g = PhysicalGrid(xlimits,ylimits,Δx)
    NX, NY = size(g)

    Fo = α*Δt/Δx^2

    Sc = Nodes{Primal,NX,NY,Float64}()
    Vf = Edges{Primal,NX,NY,Float64}()
    gradTf = Edges{Primal,NX,NY,Float64}()


    #=
    # Set up buffers
    Vf = Edges{Primal,NX,NY,Float64}()
    Vv = Edges{Primal,NX,NY,Float64}()
    Vn = Edges{Primal,NX,NY,Float64}()
    Wn = Nodes{Dual,NX,NY,Float64}()
    Vtf = EdgeGradient{Primal,Dual,NX,NY,Float64}()
    DVf = EdgeGradient{Primal,Dual,NX,NY,Float64}()
    VDVf = EdgeGradient{Primal,Dual,NX,NY,Float64}()
    =#

    #L = plan_laplacian(Sn,with_inverse=true)

    qforcefields = _process_forcing(qforce,Sc,g)


    # for now, if there are any bodies that are Open,
    # then force problem_side to ExternalInternalProblem.
    # but should be more flexible here
    problem_side_internal = _any_open_bodies(bodies) ? ExternalInternalProblem : problem_side

    N = numpts(bodies)

    Sb = ScalarData(N)
    Vb = VectorData(N)
    gradTb = VectorData(N)
    ΔTs = ScalarData(N)
    σ = ScalarData(N)
    τ = ScalarData(N)
    Cc = nothing

     points, normals, dlc, Rc, Ec, Rf, Ef, Cc =
            _immersion_operators(bodies,g,problem_side_internal,bctype,ddftype,Sc,Sb,Vf,Vb)

    conduction_L = plan_laplacian(Sc,factor=α/Δx^2)

    if isnothing(bodies)
      state_prototype = solvector(state=Sc)
      f = ConstrainedODEFunction(heatconduction_rhs!,conduction_L,_func_cache=state_prototype)
    else
      if static_points
        if bctype <: DirichletBC
          state_prototype = solvector(state=Sc,constraint=σ)
          f = ConstrainedODEFunction(heatconduction_rhs!,temperature_bc_constraint_rhs!,
                                      temperature_constraint_force!,temperature_constraint_op!,
                                      conduction_L,_func_cache=state_prototype)
        elseif bctype <: NeumannBC
          state_prototype = solvector(state=Sc)
          f = ConstrainedODEFunction(heatconduction_rhs!,conduction_L,_func_cache=state_prototype)
        end
      #=
        else
        state_prototype = solvector(state=Sc,constraint=σ,aux_state=zero_body_state(bodies))
        rhs! = ConstrainedSystems.r1vector(state_r1 = ns_rhs!,aux_r1 = rigid_body_rhs!)
        f = ConstrainedODEFunction(rhs!,temperature_bc_constraint_rhs!,
                                      temperature_constraint_force!,temperature_constraint_op!,
                                      conduction_L,_func_cache=state_prototype,
                                      param_update_func=update_immersion_operators!)
        =#
      end

    end



    HeatConduction{NX, NY, N, _motiontype(static_points), problem_side_internal, ddftype, typeof(f), typeof(state_prototype)}(
                          params, bodies, bctype, bodytemps, qforce,
                          _heat_flux_regions(qflux,Sc,g),
                          _heat_model_regions(qmodel,Sc,g),
                          _line_source(qline,Sc,g),
                          g, Δt, # rk,
                          dlc,
                          points, normals, Rc, Ec, Rf, Ef, Cc,
                          f,state_prototype,Sc,Sb,Vf,gradTf,Vb,gradTb,ΔTs,σ,τ)

end

#=
HeatConduction(Re,Δx,xlim,ylim,Δt,bodies::BodyList;
        motions=RigidMotionList(map(x -> RigidBodyMotion(0.0,0.0),bodies)),kwargs...) =
        HeatConduction(Re,Δx,xlim,ylim,Δt;bodies=bodies,motions=motions,kwargs...)

HeatConduction(Re,Δx,xlim,ylim,Δt,body::Body;kwargs...) =
        HeatConduction(Re,Δx,xlim,ylim,Δt,BodyList([body]);kwargs...)

function HeatConduction(Re,Δx,xlim,ylim,Δt,bodies::BodyList,motions::RigidMotionList;static_points=false,kwargs...)
    length(bodies) == length(motions) || error("Inconsistent lengths of bodies and motions lists")
    HeatConduction(Re,Δx,xlim,ylim,Δt,bodies;motions=motions,static_points=static_points,kwargs...)
end

HeatConduction(Re,Δx,xlim,ylim,Δt,body::Body,motion::RigidBodyMotion;static_points=false,kwargs...) =
        HeatConduction(Re,Δx,xlim,ylim,Δt,BodyList([body]),RigidMotionList([motion]);static_points=static_points,kwargs...)

=#

function Base.show(io::IO, sys::HeatConduction{NX,NY,N,MT,SD}) where {NX,NY,N,MT,SD}
    mtype = (MT == StaticPoints) ? "static" : "moving"
    sdmsg = (N == 0) ? "Unbounded" : ((SD == ExternalProblem) ? "External problem" : ((SD == InternalProblem) ? "Internal problem" : "External/internal"))
    println(io, "$sdmsg Heat conduction system on a grid of size $NX x $NY and $N $mtype immersed points")
    if N > 0
      bdmsg = (length(sys.bodies) == 1) ? "1 body" : "$(length(sys.bodies)) bodies"
      println(io, "   $bdmsg")
    end
end

# Routines to set up the immersion operators

function _immersion_operators(bodies::BodyList,g::PhysicalGrid,problem_side::Type{SD},bctype::Type{BC},ddftype,
                                Sc::ScalarGridData{NX,NY},Sb::ScalarData{N},Vf::VectorGridData{NX,NY},Vb::VectorData{N}) where {NX,NY,N,SD<:ComputationalHeatTransfer.ProblemSide,BC<:AbstractBC}

  points = VectorData(collect(bodies))
  numpts(bodies) == N || error("Inconsistent size of bodies")

  body_areas = areas(bodies)
  body_normals = normals(bodies)

  if !(problem_side==ExternalInternalProblem)
    dlc = DoubleLayer(bodies,g,Sc)
  else
    dlc = nothing
  end

  regop = _regularization(points,g,bodies,ddftype)

  Rc = RegularizationMatrix(regop,Sb,Sc) # Used by B₁ᵀ in temperature constraints
  Ec = InterpolationMatrix(regop,Sc,Sb) # Used by constraint_rhs! and B₂
  Rf = RegularizationMatrix(regop,Vb,Vf) # Used by B₁ᵀ in heat flux constraints
  Ef = InterpolationMatrix(regop,Vf,Vb) # Used by constraint_rhs! and B₂

  Cc = _neumann_matrix(bctype,body_normals,Rf,Ef,Sc,Sb,Vf,Vb)

  return points, body_normals, dlc, Rc, Ec, Rf, Ef, Cc
end


_immersion_operators(::Nothing,a...) =
    VectorData(0), VectorData(0), nothing, nothing, nothing, nothing, nothing, nothing


# For updating the system with body data

function update_immersion_operators!(sys::HeatConduction{NX,NY,N,MT,SD,DDF},bodies::BodyList) where {NX,NY,N,MT,SD,DDF<:CartesianGrids.DDFType}
    sys.bodies = deepcopy(bodies)
    sys.points, sys.dlc, sys.Rc, sys.Ec, sys.Rf, sys.Ef =
      _immersion_operators(sys.bodies,sys.grid,SD,DDF,sys.Sc,sys.Sb)
    return sys
end

update_immersion_operators!(sys::HeatConduction,body::Body) = update_immersion_operators!(sys,BodyList([body]))

function update_immersion_operators!(sys::HeatConduction,x::AbstractVector)
    tl! = RigidTransformList(x)
    tl!(sys.bodies)
    update_immersion_operators!(sys,sys.bodies)
end

# The form passed to ConstrainedODEFunction
update_immersion_operators!(sys::HeatConduction,u,sys_old::HeatConduction,t) =
    update_immersion_operators!(sys,aux_state(u))


"""
    setstepsizes(α[,gridPe=2][,cfl=0.5][,fourier=0.5]) -> Float64, Float64

Set the grid cell spacing and time step size based on the thermal diffusivity `α`,
the grid Peclet number `gridPe`, cfl number `cfl`, and grid Fourier number `fourier`.
The last three parameters all have default values.

# Example

Here is an example of setting parameters based on Reynolds number 100 (with
  default choices for grid Reynolds number, CFL number, and Fourier number):
```jldoctest
julia> Δx, Δt = setstepsizes(100)
(0.02, 0.01)
```
"""
function setstepsizes(α::Real; gridPe = 2.0, cfl = 0.5, fourier = 0.5)
    Δx = α*gridPe
    Δt = fourier*Δx^2
    return Δx, Δt
end




# some convenience functions
"""
    size(sys::HeatConduction,d::Int) -> Int

Return the number of indices of the grid used by `sys` along dimension `d`.
"""
size(sys::HeatConduction{NX,NY},d::Int) where {NX,NY} = d == 1 ? NX : NY

"""
    size(sys::HeatConduction) -> Tuple{Int,Int}

Return a tuple of the number of indices of the grid used by `sys`
"""
size(sys::HeatConduction{NX,NY}) where {NX,NY} = (size(sys,1),size(sys,2))

"""
    cellsize(sys::HeatConduction) -> Float64

Return the grid cell size of system `sys`
"""
cellsize(sys::HeatConduction) = cellsize(sys.grid)

"""
    timestep(sys::HeatConduction) -> Float64

Return the time step size of system `sys`
"""
timestep(sys::HeatConduction) = sys.Δt

"""
    origin(sys::HeatConduction) -> Tuple{Int,Int}

Return a tuple of the indices of the primal node that corresponds to the
physical origin of the coordinate system used by `sys`. Note that these
indices need not lie inside the range of indices occupied by the grid.
For example, if the range of physical coordinates occupied by the grid
is (1.0,3.0) x (2.0,4.0), then the origin is not inside the grid.
"""
origin(sys::HeatConduction) = origin(sys.grid)


"""
    timerange(tf,sys::HeatConduction)

Create a range of times, starting at the t = Δt (the time step of `sys`),
and ending at t = `tf`.
"""
timerange(tf,sys) = timestep(sys):timestep(sys):tf

# Wrap the output of the motion evaluation in VectorData
@inline assign_velocity!(u::VectorData,a...) = (assign_velocity!(u.u,u.v,a...); u)


@inline normals(sys::HeatConduction) = (!isnothing(sys.bodies)) ? normals(sys.bodies) : nothing


"""
    newstate(sys::HeatConduction)

Return a new (zero) instance of the state vector for `sys`.
"""
newstate(sys::HeatConduction) = zero(sys.state_prototype)

"""
    newstate(s::AbstractSpatialField,sys::HeatConduction)

Return an instance of the state vector for `sys`, assigned the
data in the spatial field `s`.
"""
function newstate(s::AbstractSpatialField,sys::HeatConduction)
  u = newstate(sys)
  gf = GeneratedField(state(u),s,sys.grid)
  state(u) .= cellsize(sys)*gf()
  return u
end

# Other functions
_motiontype(isstatic::Bool) = isstatic ? StaticPoints : MovingPoints

_body_closure_type(b::T) where {T<:Body{N,C}} where {N,C} = C

_any_open_bodies(nothing) = false
_any_open_bodies(bodies::BodyList) =  any(b -> _body_closure_type(b) == RigidBodyTools.OpenBody,bodies)

_regularization(sys::HeatConduction{NX, NY, N, MT, SD, DDF}) where {NX,NY,N,MT,SD,DDF} =
        _regularization(sys.points,sys.grid,sys.bodies,DDF)

_regularization(points,g,bodies,ddftype) = Regularize(points,cellsize(g),
                                I0=CartesianGrids.origin(g),weights=areas(bodies).data,ddftype=ddftype)


function _heat_flux_regions(qparams::Vector{<:PrescribedHeatFluxRegion},u::ScalarGridData,g::PhysicalGrid)
  qflux = PrescribedHeatFlux[]
  for qp in qparams
    push!(qflux,PrescribedHeatFlux(qp,u,g))
  end
  qflux
end

_heat_flux_regions(qparams::PrescribedHeatFluxRegion,u::ScalarGridData,g::PhysicalGrid) =
    _heat_flux_regions([qparams],u,g)

_heat_flux_regions(::Nothing,u,g) = nothing


function _heat_model_regions(qparams::Vector{<:PrescribedHeatModelRegion},u::ScalarGridData,g::PhysicalGrid)
  qmodel = PrescribedHeatModel[]
  for qp in qparams
    push!(qmodel,PrescribedHeatModel(qp,u,g))
  end
  qmodel
end

_heat_model_regions(qparams::PrescribedHeatModelRegion,u::ScalarGridData,g::PhysicalGrid) =
    _heat_model_regions([qparams],u,g)


_heat_model_regions(::Nothing,u,g) = nothing

function _line_source(lineparams::Vector{<:LineSourceParams},u::ScalarGridData,g::PhysicalGrid)
  qline = PrescribedLineSource[]
  for p in lineparams
    push!(qline,PrescribedLineSource(p,u,g))
  end
  qline
end

_line_source(lineparams::LineSourceParams,u::ScalarGridData,g::PhysicalGrid) =
    _line_source([lineparams],u,g)

_line_source(::Nothing,u,g) = nothing

function set_linesource_strength!(sys::HeatConduction,q::Vector{T}) where {T<:Real}
    @unpack qline = sys
    if isnothing(qline)
      return sys
    else
      total_len = mapreduce(ql -> length(ql.q),+,qline)
      @assert total_len == length(q)
      qcpy = copy(q)
      first_index = 0
      for ql in qline
        ql.q .= view(q,first_index+1:first_index+length(ql.q))
        first_index += length(ql.q)
      end
    end
    return sys
end

# Form the matrix n.Ef*Rf*n for implementing Neumann conditions
_neumann_matrix(a...) = nothing

function _neumann_matrix(::Type{<:NeumannBC},normals,Rf::RegularizationMatrix,Ef::InterpolationMatrix,
                          Sc::ScalarGridData{NX,NY},Sb::ScalarData{N},
                          Vf::VectorGridData{NX,NY},Vb::VectorData{N}) where {NX,NY,N}

  M = Matrix{Float64}(undef,N,N)

  fill!(Sb,0.0)
  for i in 1:N
    Sb[i] = 1.0
    product!(Vb,normals,Sb)
    Vf .= Rf*Vb
    Vb .= Ef*Vf
    pointwise_dot!(Sb,Vb,normals)
    M[:,i] .= Sb
    fill!(Sb,0.0)
  end
  return M
end


include("operators/surfacetemperatures.jl")
include("operators/fields.jl")
#include("operators/pointforce.jl")
include("operators/basicoperators.jl")
include("operators/surfaceoperators.jl")
include("operators/movingbodyoperators.jl")
include("operators/timemarching.jl")

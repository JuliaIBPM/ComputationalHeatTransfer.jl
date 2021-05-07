### Computing fields and other physical quantities ###

"""
    temperature(sol::ODESolution,sys::HeatConduction,t)

Return the vorticity field associated with solution vector `sol` on the grid in `sys`,
at time(s) `t`.
"""
@inline temperature(T::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY},t::Real) where {NX,NY} = T


function temperature_linesource(T::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY},t::Real;index=1) where {NX,NY}
    sys.qline[index].E*T
end



"""
    temperature(integrator)

Return the temperature field associated with `integrator` at its current state.
""" function temperature end

"""
    temperature_linesource(integrator,[index=1])

Using the solution associated with `integrator` at its current state, evaluate
the temperature field along line source `index` (which defaults to 1).
""" function temperature_linesource end



for fcn in (:temperature,:temperature_linesource)

  @eval $fcn(s::ConstrainedSystems.ArrayPartition,sys::HeatConduction,t;kwargs...) = $fcn(state(s),sys,t;kwargs...)

end

for fcn in (:temperature,:temperature_linesource)
  @eval $fcn(integ::ConstrainedSystems.OrdinaryDiffEq.ODEIntegrator;kwargs...) = $fcn(integ.u,integ.p,integ.t;kwargs...)

  @eval $fcn(sol::ConstrainedSystems.OrdinaryDiffEq.ODESolution,sys::HeatConduction,t;kwargs...) = $fcn(sol(t),sys,t;kwargs...)

  @eval $fcn(sol::ConstrainedSystems.OrdinaryDiffEq.ODESolution,sys::HeatConduction,t::AbstractArray;kwargs...) =
      map(ti -> $fcn(sol(ti),sys,ti;kwargs...),t)

end




# Surface field quantities

#=
"""
    traction(integrator)

Return the traction at all of the surface points associated with `integrator` at
its current state.
""" function traction end

"""
    pressurejump(integrator)

Return the pressure jump at all of the surface points associated with `integrator` at
its current state.
""" function pressurejump end


for fcn in (:traction,:pressurejump)
  @eval $fcn(s::ConstrainedSystems.ArrayPartition,a...;kwargs...) = $fcn(constraint(s),a...;kwargs...)
  @eval $fcn(integ::ConstrainedSystems.OrdinaryDiffEq.ODEIntegrator) = $fcn(integ.u,integ.p,integ.t)
  @eval $fcn(sol::ODESolution,sys::HeatConduction,t) = $fcn(sol(t),sys,t)

  @eval $fcn(sol::ODESolution,sys::HeatConduction,t::AbstractArray) =
        map(ti -> $fcn(sol(ti),sys,ti),t)

end

"""
    traction(sol::ODESolution,sys::HeatConduction,t)

Return the traction at all of the surface points associated with solution `sol` on
grid in `sys` at time(s) t.
"""
traction(τ,sys::HeatConduction{NX,NY,0,MT,FS,SD},t) where {NX,NY,MT,FS,SD} = τ

traction(τ::VectorData{N},sys::HeatConduction{NX,NY,N,MT,FS,ExternalInternalFlow},t) where {NX,NY,N,MT,FS} = τ

function traction(τ::VectorData{N},sys::HeatConduction{NX,NY,N,MT,FS,SD},t) where {NX,NY,N,MT,FS,SD}
    nrm = normals(sys.bodies)
    relative_surface_velocity!(sys.Δus,sys,t)
    pointwise_dot!(sys.Sb,nrm,sys.Δus)
    surface_velocity_jump!(sys.Δus,sys,t)
    product!(sys.Vb,sys.Δus,sys.Sb)
    return τ + sys.Vb
end

"""
    pressurejump(sol::ODESolution,sys::HeatConduction,t)

Return the pressure jump at all of the surface points associated with solution `sol` on
grid in `sys` at time(s) t.
"""
function pressurejump(τ::VectorData{N},sys::HeatConduction{NX,NY,N},t) where {NX,NY,N}

    #_hasfilter(sys) ? (τf = similar(τ); τf .= sys.Cf*force(τ,sys)) : τf = deepcopy(force(τ,sys))
    #τf = deepcopy(τ)
    #return nrm.u∘τf.u + nrm.v∘τf.v # This might need some modification

    press = zero(sys.Sb)
    nrm = normals(sys.bodies)
    pointwise_dot!(press,nrm,traction(τ,sys,t))
    return -press
end

## Total quantities

force(τ,sys::HeatConduction{NX,NY,0},t,bodyi::Int) where {NX,NY} = Vector{Float64}(), Vector{Float64}()

function force(τ::VectorData{N},sys::HeatConduction,t,bodyi::Int) where {N}
    product!(sys.τ,traction(τ,sys,t),areas(sys.bodies))
    fxvec = sum(sys.τ.u,sys.bodies,bodyi)
    fyvec = sum(sys.τ.v,sys.bodies,bodyi)
    return fxvec, fyvec
end


force(s::ConstrainedSystems.ArrayPartition,sys,t,bodyi) = force(constraint(s),sys,t,bodyi)

"""
    force(integ,bodyindex)

Given the state of the system in integrator `integ`, return the current force
on the body with index `bodyindex`.
"""
force(integ::ConstrainedSystems.OrdinaryDiffEq.ODEIntegrator,bodyi) = force(integ.u,integ.p,integ.t,bodyi)

"""
    force(sol,sys::HeatConduction,bodyindex) -> Tuple{Vector,Vector}

Given the solution history vector `sol` and the system `sys`, return the force
history on the body with index `bodyindex` as a tuple of vectors, one for
each component.
"""
function force(sol::ConstrainedSystems.OrdinaryDiffEq.ODESolution,sys::HeatConduction,bodyi::Int)
    fx = map((u,t) -> force(u,sys,t,bodyi)[1],sol.u,sol.t)
    fy = map((u,t) -> force(u,sys,t,bodyi)[2],sol.u,sol.t)
    fx, fy
end


function moment(τ::VectorData{N},bodies::BodyList,sys::HeatConduction,t,bodyi::Int;center=(0.0,0.0)) where {N}
    xc, yc = center
    x, y = collect(bodies)
    dX = VectorData(x .- xc,  y .- yc)
    product!(sys.τ,traction(τ,sys,t),areas(bodies))
    return sum(cross(dX,sys.τ),bodies,bodyi)
end


function moment(s::ConstrainedSystems.ArrayPartition,sys::HeatConduction{NX,NY,N,MovingPoints},t,bodyi;center=(0.0,0.0)) where {NX,NY,N}
    bodies = deepcopy(sys.bodies)
    tl! = RigidTransformList(aux_state(s))
    tl!(bodies)
    moment(constraint(s),bodies,sys,t,bodyi;center=center)
end

moment(s::ConstrainedSystems.ArrayPartition,sys::HeatConduction{NX,NY,N,StaticPoints},t,bodyi;center=(0.0,0.0)) where {NX,NY,N} =
        moment(constraint(s),sys.bodies,sys,t,bodyi;center=center)


"""
    moment(integ,bodyindex[,center=(0,0))

Given the state of the system in integrator `integ`, return the current moment
on the body with index `bodyindex`.  The center of the
moment can be passed as an optional tuple argument.
"""
moment(integ::ConstrainedSystems.OrdinaryDiffEq.ODEIntegrator,bodyi;center=(0.0,0.0)) =
      moment(integ.u,integ.p,integ.t,bodyi,center=center)


"""
    moment(sol,sys::HeatConduction,bodyindex[,center=(0,0)]) -> Vector

Given the solution history vector `sol` and the system `sys`, return the moment
history on the body with index `bodyindex` as a vector. The center of the
moment can be passed as an optional tuple argument.
"""
function moment(sol::ConstrainedSystems.OrdinaryDiffEq.ODESolution,sys::HeatConduction,bodyi::Int;center=(0.0,0.0))
    mom = map((u,t) -> moment(u,sys,t,bodyi,center=center),sol.u,sol.t)
end
=#

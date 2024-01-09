import CartesianGrids: dot
import ConstrainedSystems: init, solve

_norm_sq(u) = dot(u,u)
_norm_sq(u::ConstrainedSystems.ArrayPartition) = sum(_norm_sq,u.x)
state_norm(u,t) = sqrt(_norm_sq(u))

# ensure that time marching makes use of
#init(prob;alg=LiskaIFHERK(),kwargs...) = init(prob,alg;internal_norm=state_norm,kwargs...)

# function init(u0,tspan,sys::HeatConduction;alg=ConstrainedSystems.LiskaIFHERK(),kwargs...)
#     prob = ODEProblem(sys.f,u0,tspan,sys)
#     return init(prob, alg,dt=timestep(sys),internal_norm=state_norm,kwargs...)
# end

# function solve(prob,sys::HeatConduction;alg=ConstrainedSystems.LiskaIFHERK(),kwargs...)
#   return solve(prob,alg,dt=timestep(sys),internal_norm=state_norm,kwargs...)
# end

function init(u0,tspan,sys::HeatConduction;alg=ConstrainedSystems.LiskaIFHERK(),kwargs...)
  prob = ODEProblem(sys.f,u0,tspan,sys)
  return init(prob, alg,dt=timestep(sys),kwargs...)
end

function solve(prob,sys::HeatConduction;alg=ConstrainedSystems.LiskaIFHERK(),kwargs...)
return solve(prob,alg,dt=timestep(sys),kwargs...)
end

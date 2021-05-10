module ComputationalHeatTransfer

using DocStringExtensions
using Reexport

@reexport using CartesianGrids
@reexport using ImmersedLayers
@reexport using RigidBodyTools
@reexport using ConstrainedSystems
@reexport using GridUtilities


using LinearAlgebra
#using SparseArrays

export HeatConduction,setstepsizes,timestep,newstate,
       ExternalProblem, InternalProblem, ExternalInternalProblem,
       temperature,temperature_linesource,
       HeatConductionParameters,LineSourceParams,
       PrescribedHeatFluxRegion,PrescribedHeatModelRegion,
       set_linesource_strength!

#=
export NavierStokes, PulseParams, PointForce, SpatialDGaussian,
       setstepsizes, timestep, timerange, newstate,
       update_immersion_operators!,
       vorticity, velocity, velocity!, streamfunction, streamfunction!,
       scalarpotential, scalarpotential!, convective_derivative, convective_derivative!,
       pressure, traction, force, moment, pressurejump
=#

abstract type PointMotionType end
abstract type StaticPoints <: PointMotionType end
abstract type MovingPoints <: PointMotionType end

abstract type ProblemSide end
abstract type ExternalProblem <: ProblemSide end
abstract type InternalProblem <: ProblemSide end
abstract type ExternalInternalProblem <: ProblemSide end

const NDIM = 2

include("utils/forcing.jl")
include("utils/ohp.jl")
include("heat_conduction.jl")
#include("plot_recipes.jl")



end

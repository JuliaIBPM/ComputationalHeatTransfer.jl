module ComputationalHeatTransfer

using DocStringExtensions
using Reexport
# using Revise
using LinearAlgebra

@reexport using CartesianGrids
@reexport using ImmersedLayers
@reexport using RigidBodyTools
@reexport using ConstrainedSystems
@reexport using GridUtilities

export HeatConduction,setstepsizes,timestep,newstate,
       ExternalProblem, InternalProblem, ExternalInternalProblem,
       temperature,temperature_linesource,
       HeatConductionParameters,LineSourceParams,
       PrescribedHeatFluxRegion,PrescribedHeatModelRegion,
       set_linesource_strength!



abstract type PointMotionType end
abstract type StaticPoints <: PointMotionType end
abstract type MovingPoints <: PointMotionType end

abstract type ProblemSide end
abstract type ExternalProblem <: ProblemSide end
abstract type InternalProblem <: ProblemSide end
abstract type ExternalInternalProblem <: ProblemSide end

abstract type AbstractBC end
abstract type NeumannBC <: AbstractBC end
abstract type DirichletBC <: AbstractBC end
abstract type AdiabaticBC <: NeumannBC end
abstract type HomogeneousDirichletBC <: DirichletBC end


const NDIM = 2

include("utils/forcing.jl")
# include("utils/ohp.jl")
include("heat_conduction.jl")
# include("OneDOHP.jl")
# include("fileIO.jl")

end

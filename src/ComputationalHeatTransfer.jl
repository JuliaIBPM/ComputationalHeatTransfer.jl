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



include("defaults.jl")
include("fields.jl")
include("ode_operators.jl")
include("api.jl")

end

#####################################################################


# ComputationalHeatTransfer.jl
tools for numerical simulation of conductive and convective heat transfer

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaIBPM.github.io/ComputationalHeatTransfer.jl/stable) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaIBPM.github.io/ComputationalHeatTransfer.jl/dev)
[![Build Status](https://github.com/JuliaIBPM/ComputationalHeatTransfer.jl/workflows/CI/badge.svg)](https://github.com/JuliaIBPM/ComputationalHeatTransfer.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaIBPM/ComputationalHeatTransfer.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaIBPM/ComputationalHeatTransfer.jl)



## About the package

The purpose of this package is to enable easy setup and solution of steady and transient heat conduction and convection problems in arbitrary geometries.  Documentation can be found at https://JuliaIBPM.github.io/ComputationalHeatTransfer.jl/latest.

**ComputationalHeatTransfer.jl** is registered in the general Julia registry. To install, type
e.g.,
```julia
] add ComputationalHeatTransfer
```

Then, in any version, type
```julia
julia> using ComputationalHeatTransfer
```
For examples, consult the documentation or see the example Jupyter notebooks in the Examples folder.

Many of the core aspects of the fluid-body interaction are based on the Method of Immersed Layers [1], which is an extension of the immersed boundary projection method [2].

[1]: Eldredge, J. D. (2022) "A method of immersed layers on Cartesian grids, with application to incompressible flows," *J. Comput. Phys.*, 448, [110716](https://arxiv.org/abs/2103.04521).


[2]: Taira, K. and Colonius, T. (2007) "The immersed boundary method: a projection approach," *J. Comput. Phys.*, 225, 2118--2137.

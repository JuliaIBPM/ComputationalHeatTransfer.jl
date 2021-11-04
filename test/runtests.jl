# Pkg.clone("https://github.com/CoolProp/CoolProp.jl.git")
# Pkg.build("CoolProp")

# using ComputationalHeatTransfer
# using LaTeXStrings
using Revise
# using LinearAlgebra
# using DifferentialEquations
# using Interpolations
# using JLD2


cd(dirname(pwd()))
cd("src")
includet("OneDOHP.jl")
using ..OneDOHP

using Test

@testset "trigonometric identities" begin
           θ = 2/3*π
           @test sin(-θ) ≈ -sin(θ)
           @test cos(-θ) ≈ cos(θ)
           @test sin(2θ) ≈ 2*sin(θ)*cos(θ)
           @test cos(2θ) ≈ cos(θ)^2 - sin(θ)^2
       end;

# ##using TestSetExtensions
# using Literate
#
# const GROUP = get(ENV, "GROUP", "All")
#
# notebookdir = "../examples"
# docdir = "../docs/src/manual"
# litdir = "./literate"
#
# if GROUP == "All" || GROUP == "Auxiliary"
#   #include("pointforce.jl")
# end
#
#
# if GROUP == "All" || GROUP == "Notebooks"
#   for (root, dirs, files) in walkdir(litdir)
#     for file in files
#       #endswith(file,".jl") && startswith(file,"6") && Literate.notebook(joinpath(root, file),notebookdir)
#       endswith(file,".jl") && Literate.notebook(joinpath(root, file),notebookdir)
#     end
#   end
# end
#
# if GROUP == "Documentation"
#   for (root, dirs, files) in walkdir(litdir)
#     for file in files
#       endswith(file,".jl") && Literate.markdown(joinpath(root, file),docdir)
#     end
#   end
# end

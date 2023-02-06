using Pkg
Pkg.activate(dirname(pwd()))

using Documenter
using ComputationalHeatTransfer

makedocs(
    prettyurls = true,
    sitename = "ComputationalHeatTransfer",
    format = Documenter.HTML(),
    modules = [ComputationalHeatTransfer],
    pages = [
        "Home" => "index.md",
        "Basics" => ["manual/subpage.md"
                     ]
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "<github.com/JuliaIBPM/ComputationalHeatTransfer.jl.git>",
    target = "build",
    deps = nothing,
    make = nothing,
    versions = nothing
)


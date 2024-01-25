using Documenter, ViscousFlow

ENV["GKSwstype"] = "nul" # removes GKS warnings during plotting

makedocs(
    sitename = "ViscousFlow.jl",
    doctest = true,
    clean = true,
    modules = [ViscousFlow],
    pages = [
        "Home" => "index.md",
        "Manual" => ["manual/heatconduction_dirichlet.md",
                    "manual/heatconduction_neumann.md",
                    "manual/heatconduction_unbounded.md"
                     ]
        #"Internals" => [ "internals/properties.md"]
    ],
    #format = Documenter.HTML(assets = ["assets/custom.css"])
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine = MathJax(Dict(
            :TeX => Dict(
                :equationNumbers => Dict(:autoNumber => "AMS"),
                :Macros => Dict()
            )
        ))
    ),
    #assets = ["assets/custom.css"],
    #strict = true
)


#if "DOCUMENTER_KEY" in keys(ENV)
deploydocs(
     repo = "github.com/JuliaIBPM/ViscousFlow.jl.git",
     target = "build",
     deps = nothing,
     make = nothing
     #versions = "v^"
)
#end

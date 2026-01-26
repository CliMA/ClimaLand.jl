using SIMDTypes
using Documenter

DocMeta.setdocmeta!(SIMDTypes, :DocTestSetup, :(using SIMDTypes); recursive=true)

makedocs(;
    modules=[SIMDTypes],
    authors="Julia Computing, Inc. and Contributors",
    repo="https://github.com/JuliaSIMD/SIMDTypes.jl/blob/{commit}{path}#{line}",
    sitename="SIMDTypes.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaSIMD.github.io/SIMDTypes.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaSIMD/SIMDTypes.jl",
)

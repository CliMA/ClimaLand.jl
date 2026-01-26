using BitTwiddlingConvenienceFunctions
using Documenter

DocMeta.setdocmeta!(BitTwiddlingConvenienceFunctions, :DocTestSetup, :(using BitTwiddlingConvenienceFunctions); recursive=true)

makedocs(;
    modules=[BitTwiddlingConvenienceFunctions],
    authors="Julia Computing",
    repo="https://github.com/JuliaSIMD/BitTwiddlingConvenienceFunctions.jl/blob/{commit}{path}#{line}",
    sitename="BitTwiddlingConvenienceFunctions.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaSIMD.github.io/BitTwiddlingConvenienceFunctions.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaSIMD/BitTwiddlingConvenienceFunctions.jl",
    devbranch="main",
)

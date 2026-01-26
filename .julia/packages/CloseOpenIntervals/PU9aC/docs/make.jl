using CloseOpenIntervals
using Documenter

DocMeta.setdocmeta!(CloseOpenIntervals, :DocTestSetup, :(using CloseOpenIntervals); recursive=true)

makedocs(;
    modules=[CloseOpenIntervals],
    authors="chriselrod <elrodc@gmail.com> and contributors",
    repo="https://github.com/JuliaSIMD/CloseOpenIntervals.jl/blob/{commit}{path}#{line}",
    sitename="CloseOpenIntervals.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaSIMD.github.io/CloseOpenIntervals.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaSIMD/CloseOpenIntervals.jl",
    devbranch = "main",
)

using FunctionWrappersWrappers
using Documenter

DocMeta.setdocmeta!(FunctionWrappersWrappers, :DocTestSetup, :(using FunctionWrappersWrappers); recursive=true)

makedocs(;
    modules=[FunctionWrappersWrappers],
    authors="Chris Elrod <elrodc@gmail.com> and contributors",
    repo="https://github.com/chriselrod/FunctionWrappersWrappers.jl/blob/{commit}{path}#{line}",
    sitename="FunctionWrappersWrappers.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://chriselrod.github.io/FunctionWrappersWrappers.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/chriselrod/FunctionWrappersWrappers.jl",
    devbranch="main",
)

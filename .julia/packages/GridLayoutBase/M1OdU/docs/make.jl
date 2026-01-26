using Documenter, GridLayoutBase

makedocs(;
    modules=[GridLayoutBase],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/jkrumbiegel/GridLayoutBase.jl/blob/{commit}{path}#L{line}",
    sitename="GridLayoutBase.jl",
    authors="Julius Krumbiegel",
    assets=String[],
)

deploydocs(;
    repo="github.com/jkrumbiegel/GridLayoutBase.jl",
)

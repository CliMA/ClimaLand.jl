using Documenter, ConstructionBase

makedocs(;
    modules=[ConstructionBase],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/JuliaObjects/ConstructionBase.jl/blob/{commit}{path}#L{line}",
    sitename="ConstructionBase.jl",
    authors="Takafumi Arakaki, Rafael Schouten, Jan Weidner",
    strict=true,
)

deploydocs(;
    repo="github.com/JuliaObjects/ConstructionBase.jl",
)

using Documenter, PolygonOps

makedocs(;
    modules=[PolygonOps],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/JuliaGeometry/PolygonOps.jl/blob/{commit}{path}#L{line}",
    sitename="PolygonOps.jl",
    authors="steve <kd2cca@gmail.com>")

deploydocs(;
    repo="github.com/JuliaGeometry/PolygonOps.jl",
)

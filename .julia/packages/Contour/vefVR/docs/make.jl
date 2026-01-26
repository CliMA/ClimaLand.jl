using Documenter, Contour

makedocs(
    modules = [Contour],
    format = Documenter.HTML(),
    sitename = "Contour.jl",
    pages = Any[
        "Introduction" => "index.md",
        "Tutorial" => "tutorial.md",
        "Reference" => "reference.md",
    ]
)

deploydocs(
    repo = "github.com/JuliaGeometry/Contour.jl.git",
)

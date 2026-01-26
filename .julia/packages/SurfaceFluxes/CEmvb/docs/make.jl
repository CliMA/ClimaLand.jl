using SurfaceFluxes, Documenter
using DocumenterCitations

# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "nul"

bib = CitationBibliography(joinpath(@__DIR__, "bibliography.bib"))

#! format: off
pages = Any[
    "Home" => "index.md",
    "References" => "References.md",
    "Equations" => "SurfaceFluxes.md",
    "Universal Functions" => "UniversalFunctions.md"
]

mathengine = MathJax(Dict(
    :TeX => Dict(
        :equationNumbers => Dict(:autoNumber => "AMS"),
        :Macros => Dict(),
    ),
))

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    mathengine = mathengine,
    collapselevel = 1,
)
#! format: on

makedocs(
    bib,
    sitename = "SurfaceFluxes.jl",
    strict = true,
    format = format,
    checkdocs = :exports,
    clean = true,
    doctest = true,
    modules = [SurfaceFluxes],
    pages = pages,
)

deploydocs(
    repo = "github.com/CliMA/SurfaceFluxes.jl.git",
    target = "build",
    push_preview = true,
    devbranch = "main",
    forcepush = true,
)

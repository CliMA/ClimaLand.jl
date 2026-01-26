using Documenter
using ClimaDiagnostics
import ClimaUtilities

pages = [
    "Overview" => "index.md",
    "User guide" => "user_guide.md",
    "Saving output" => "writers.md",
    "How to add ClimaDiagnostics to a package" => "developer_guide.md",
    "Internals" => "internals.md",
    "APIs" => "api.md",
]

mathengine = MathJax(
    Dict(
        :TeX => Dict(
            :equationNumbers => Dict(:autoNumber => "AMS"),
            :Macros => Dict(),
        ),
    ),
)
format = Documenter.HTML(
    prettyurls = !isempty(get(ENV, "CI", "")),
    collapselevel = 1,
    mathengine = mathengine,
)

DocMeta.setdocmeta!(
    ClimaDiagnostics,
    :DocTestSetup,
    :(using Dates);
    recursive = true,
)

makedocs(
    sitename = "ClimaDiagnostics.jl",
    authors = "Gabriele Bozzola",
    format = format,
    pages = pages,
    checkdocs = :exports,
    doctest = true,
    strict = false,
    clean = true,
    modules = [ClimaDiagnostics],
)

deploydocs(
    repo = "github.com/CliMA/ClimaDiagnostics.jl.git",
    target = "build",
    push_preview = true,
    devbranch = "main",
    forcepush = true,
)

push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using Documenter
using ClimaLSM

pages = Any["Home" => "index.md", "Contribution guide" => "Contributing.md"]

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

makedocs(
    sitename = "ClimaLSM.jl",
    authors = "Clima Land Model Team",
    format = format,
    pages = pages,
    checkdocs = :exports,
    doctest = true,
    strict = false,
    clean = true,
    modules = [ClimaLSM],
)

deploydocs(
    repo = "github.com/CliMA/ClimaLSM.jl.git",
    target = "build",
    push_preview = true,
    devbranch = "main",
    forcepush = true,
)

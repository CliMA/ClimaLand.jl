using ClimaParams, Documenter

pages = Any[
    "Home" => "index.md",
    "TOML file interface" => "toml.md",
    "Parameter retrieval" => "param_retrieval.md",
    "API" => "API.md",
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
    prettyurls = get(ENV, "CI", nothing) == "true",
    mathengine = mathengine,
    collapselevel = 1,
)

makedocs(
    sitename = "ClimaParams.jl",
    format = format,
    clean = true,
    checkdocs = :exports,
    strict = true,
    modules = [ClimaParams],
    pages = pages,
)

deploydocs(
    repo = "github.com/CliMA/ClimaParams.jl.git",
    target = "build",
    push_preview = true,
    devbranch = "main",
    forcepush = true,
)

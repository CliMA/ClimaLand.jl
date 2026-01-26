using RootSolvers, Documenter

pages = Any[
    "Home" => "index.md",
    "Getting Started" => "GettingStarted.md",
    "API" => "API.md",
    "Developer Documentation" => "DeveloperDocs.md",
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

makedocs(;
    sitename = "RootSolvers.jl",
    format = format,
    clean = true,
    modules = [RootSolvers],
    pages = pages,
)

deploydocs(
    repo = "github.com/CliMA/RootSolvers.jl.git",
    target = "build",
    push_preview = true,
)

push!(LOAD_PATH, "../src/")

using Documenter, ExactPredicates

makedocs(
    sitename="ExactPredicates.jl",
    pages = ["index.md", "api.md"],
    repo = "https://github.com/lairez/ExactPredicates.jl/blob/master{path}#{line}"
)

deploydocs(
    repo = "github.com/lairez/ExactPredicates.jl.git",
    devbranch = "master",
    versions = ["stable" => "v^"]
)

using Documenter, RoundingEmulator

makedocs(
    sitename = "RoundingEmulator.jl",
    pages = [
        "Home" => "index.md",
        "Functions" => "functions.md",
        "References" => "references.md"
    ]
)

deploydocs(repo = "github.com/matsueushi/RoundingEmulator.jl")

using Documenter, DiffResults

makedocs(
    modules=[DiffResults],
    doctest = false,
    sitename = "DiffResults",
    pages = ["Documentation" => "index.md"])

deploydocs(
    repo = "github.com/JuliaDiff/DiffResults.jl.git")

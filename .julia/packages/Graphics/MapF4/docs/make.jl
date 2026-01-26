using Documenter
using Graphics

DocMeta.setdocmeta!(Graphics, :DocTestSetup, :(using Graphics); recursive=true)

makedocs(
    sitename = "Graphics",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    modules = [Graphics],
    pages = ["index.md", "reference.md"],
    clean = false,
    checkdocs = :exports
)

deploydocs(
    repo = "github.com/JuliaGraphics/Graphics.jl.git",
    push_preview = true
)

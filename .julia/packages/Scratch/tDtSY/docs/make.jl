using Documenter, Scratch

makedocs(
    modules = [Scratch],
    sitename = "Scratch.jl",
)

deploydocs(
    repo = "github.com/JuliaPackaging/Scratch.jl.git",
    push_preview = true,
)

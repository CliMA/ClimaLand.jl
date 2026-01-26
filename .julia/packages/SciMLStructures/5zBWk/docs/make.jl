using Documenter, SciMLStructures

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

pages = [
    "Home" => "index.md",
    "interface.md",
    "example.md",
    "api.md"
]

ENV["GKSwstype"] = "100"

makedocs(modules = [SciMLStructures],
    sitename = "SciMLStructures.jl",
    clean = true,
    doctest = false,
    linkcheck = true,
    format = Documenter.HTML(assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/SciMLStructures/stable/"),
    pages = pages)

deploydocs(repo = "github.com/SciML/SciMLStructures.jl"; push_preview = true)

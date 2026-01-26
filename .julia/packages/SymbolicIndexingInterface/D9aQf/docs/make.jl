using Documenter, SymbolicIndexingInterface

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

include("pages.jl")

makedocs(sitename = "SymbolicIndexingInterface.jl",
    authors = "Chris Rackauckas",
    modules = [SymbolicIndexingInterface],
    clean = true, doctest = false,
    format = Documenter.HTML(analytics = "UA-90474609-3",
        assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/SymbolicIndexingInterface/stable/"),
    pages = pages,
    checkdocs = :exports)

deploydocs(repo = "github.com/SciML/SymbolicIndexingInterface.jl.git";
    push_preview = true)

using Documenter, SciMLOperators

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

include("pages.jl")

makedocs(sitename = "SciMLOperators.jl",
    authors = "Vedant Puri, Alex Jones, Chris Rackauckas",
    modules = [SciMLOperators],
    clean = true, doctest = false, linkcheck = true,
    warnonly = [:docs_block, :missing_docs, :cross_references, :linkcheck],
    format = Documenter.HTML(assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/SciMLOperators/stable"),
    pages = pages)

deploydocs(repo = "github.com/SciML/SciMLOperators.jl.git";
    push_preview = true)

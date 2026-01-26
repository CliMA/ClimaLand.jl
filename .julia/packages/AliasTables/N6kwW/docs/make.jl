using AliasTables
using Documenter
using DocumenterVitepress

DocMeta.setdocmeta!(AliasTables, :DocTestSetup, :(using AliasTables); recursive=true)

makedocs(;
    authors="Lilith Orion Hafner <lilithhafner@gmail.com> and contributors",
    repo=Remotes.GitHub("LilithHafner", "AliasTables.jl"),
    sitename="AliasTables.jl",
    format=DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/LilithHafner/AliasTables.jl",
        devbranch = "main",
        devurl = "dev",
        deploy_url = "aliastables.lilithhafner.com"),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/LilithHafner/AliasTables.jl",
    push_preview=true,
    devbranch="main",
)

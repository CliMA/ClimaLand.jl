using Documenter
using Jieko
using DocThemeIndigo

indigo = DocThemeIndigo.install(Jieko)

makedocs(;
    modules = [Jieko],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical="https://Roger-luo.github.io/Jieko.jl",
        assets=String[indigo],
    ),
    pages = [
        "Home" => "index.md",
    ],
    repo = "https://github.com/Roger-luo/Jieko.jl",
    sitename = "Jieko.jl",
)

deploydocs(; repo = "https://github.com/Roger-luo/Jieko.jl")

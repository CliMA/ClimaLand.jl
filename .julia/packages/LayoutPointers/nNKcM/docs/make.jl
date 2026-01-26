using LayoutPointers
using Documenter

DocMeta.setdocmeta!(
  LayoutPointers,
  :DocTestSetup,
  :(using LayoutPointers);
  recursive = true,
)

makedocs(;
  modules = [LayoutPointers],
  authors = "JuliaSIMD <elrodc@gmail.com> and contributors",
  repo = "https://github.com/JuliaSIMD/LayoutPointers.jl/blob/{commit}{path}#{line}",
  sitename = "LayoutPointers.jl",
  format = Documenter.HTML(;
    prettyurls = get(ENV, "CI", "false") == "true",
    canonical = "https://JuliaSIMD.github.io/LayoutPointers.jl",
    assets = String[],
  ),
  pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/JuliaSIMD/LayoutPointers.jl")

using CPUSummary
using Documenter

DocMeta.setdocmeta!(CPUSummary, :DocTestSetup, :(using CPUSummary); recursive = true)

makedocs(;
  modules = [CPUSummary],
  authors = "Chris Elrod <elrodc@gmail.com> and contributors",
  repo = "https://github.com/JuliaSIMD/CPUSummary.jl/blob/{commit}{path}#{line}",
  sitename = "CPUSummary.jl",
  format = Documenter.HTML(;
    prettyurls = get(ENV, "CI", "false") == "true",
    canonical = "https://JuliaSIMD.github.io/CPUSummary.jl",
    assets = String[],
  ),
  pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/JuliaSIMD/CPUSummary.jl")

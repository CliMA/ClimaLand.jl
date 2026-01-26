using StrideArraysCore
using Documenter

makedocs(;
  modules = [StrideArraysCore],
  authors = "Chris Elrod <elrodc@gmail.com> and contributors",
  repo = "https://github.com/JuliaSIMD/StrideArraysCore.jl/blob/{commit}{path}#L{line}",
  sitename = "StrideArraysCore.jl",
  format = Documenter.HTML(;
    prettyurls = get(ENV, "CI", "false") == "true",
    canonical = "https://JuliaSIMD.github.io/StrideArraysCore.jl",
    assets = String[]
  ),
  pages = ["Home" => "index.md"]
)

deploydocs(; repo = "github.com/JuliaSIMD/StrideArraysCore.jl")

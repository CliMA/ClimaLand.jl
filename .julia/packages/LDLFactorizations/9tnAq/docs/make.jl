using Documenter, LDLFactorizations

makedocs(
  modules = [LDLFactorizations],
  doctest = true,
  linkcheck = true,
  format = Documenter.HTML(
    assets = ["assets/style.css"],
    prettyurls = get(ENV, "CI", nothing) == "true",
  ),
  sitename = "LDLFactorizations.jl",
  pages = Any["Home" => "index.md", "Tutorial" => "tutorial.md", "Reference" => "reference.md"],
)

deploydocs(
  repo = "github.com/JuliaSmoothOptimizers/LDLFactorizations.jl.git",
  push_preview = true,
  devbranch = "main",
)

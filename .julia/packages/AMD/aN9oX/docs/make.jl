using Documenter, AMD

makedocs(
  modules = [AMD],
  doctest = true,
  linkcheck = true,
  format = Documenter.HTML(
    assets = ["assets/style.css"],
    prettyurls = get(ENV, "CI", nothing) == "true",
  ),
  sitename = "AMD.jl",
  pages = Any["Home" => "index.md", "Tutorial" => "tutorial.md", "Reference" => "reference.md"],
)

deploydocs(
  repo = "github.com/JuliaSmoothOptimizers/AMD.jl.git",
  push_preview = true,
  devbranch = "main",
)

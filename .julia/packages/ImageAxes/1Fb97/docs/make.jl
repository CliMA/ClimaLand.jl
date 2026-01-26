using Documenter, ImageAxes, SimpleTraits

makedocs(modules  = [ImageAxes],
         format   = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
         sitename = "ImageAxes",
         linkcheck = !("skiplinks" in ARGS),
         pages    = ["index.md", "reference.md"])

deploydocs(
           repo   = "github.com/JuliaImages/ImageAxes.jl.git",
           )

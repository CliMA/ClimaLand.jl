using Documenter, ImageBase, ImageShow

format = Documenter.HTML(edit_link = "master",
                         prettyurls = get(ENV, "CI", nothing) == "true")

makedocs(modules  = [ImageBase, ],
         format   = format,
         sitename = "ImageBase",
         pages    = ["index.md", "reference.md"])

deploydocs(repo   = "github.com/JuliaImages/ImageBase.jl.git")

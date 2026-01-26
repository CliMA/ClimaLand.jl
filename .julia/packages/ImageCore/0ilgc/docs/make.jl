using Documenter, ImageCore, ImageShow

format = Documenter.HTML(edit_link = "master",
                         prettyurls = get(ENV, "CI", nothing) == "true")

makedocs(modules  = [ImageCore],
         format   = format,
         sitename = "ImageCore",
         pages    = ["index.md", "views.md", "map.md", "traits.md", "reference.md"],
         warnonly = [:missing_docs, :cross_references])

deploydocs(repo   = "github.com/JuliaImages/ImageCore.jl.git")

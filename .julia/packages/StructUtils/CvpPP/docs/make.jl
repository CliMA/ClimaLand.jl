using Documenter, StructUtils

makedocs(modules = [StructUtils], sitename = "StructUtils.jl")

deploydocs(repo = "github.com/JuliaServices/StructUtils.jl.git", push_preview = true)

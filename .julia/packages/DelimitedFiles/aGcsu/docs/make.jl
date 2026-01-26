using DelimitedFiles
using Documenter: DocMeta, makedocs, deploydocs

DocMeta.setdocmeta!(DelimitedFiles, :DocTestSetup, :(using DelimitedFiles); recursive=true)

makedocs(
    modules = [DelimitedFiles],
    sitename = "DelimitedFiles",
    pages = Any[
        "DelimitedFiles" => "index.md"
        ];
    # strict = true,
    strict = Symbol[:doctest],
    )

deploydocs(repo = "github.com/JuliaData/DelimitedFiles.jl.git")

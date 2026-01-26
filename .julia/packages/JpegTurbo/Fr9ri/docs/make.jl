using JpegTurbo
using Documenter

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true"
)

makedocs(;
    modules=[JpegTurbo],
    sitename="JpegTurbo.jl",
    format=format,
    pages=[
        "Home" => "index.md",
        "Reference" => "reference.md",
    ],
    doctest=false, # doctest are run in unit test
)

deploydocs(;
    repo="github.com/JuliaIO/JpegTurbo.jl",
    devbranch="master",
)

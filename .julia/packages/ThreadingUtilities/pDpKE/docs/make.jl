using ThreadingUtilities
using Documenter

makedocs(;
    modules=[ThreadingUtilities],
    authors="Chris Elrod <elrodc@gmail.com> and contributors",
    repo="https://github.com/JuliaSIMD/ThreadingUtilities.jl/blob/{commit}{path}#L{line}",
    sitename="ThreadingUtilities.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaSIMD.github.io/ThreadingUtilities.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Public API" => "public-api.md",
        "Internals (Private)" => "internals.md",
    ]
)

deploydocs(;
    repo="github.com/JuliaSIMD/ThreadingUtilities.jl",
)

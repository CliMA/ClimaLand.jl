using Isoband
using Documenter

makedocs(;
    modules=[Isoband],
    authors="Julius Krumbiegel <julius.krumbiegel@gmail.com> and contributors",
    repo="https://github.com/jkrumbiegel/Isoband.jl/blob/{commit}{path}#L{line}",
    sitename="Isoband.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jkrumbiegel.github.io/Isoband.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jkrumbiegel/Isoband.jl",
)

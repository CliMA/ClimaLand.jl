using Documenter
using Atomix

makedocs(;
    sitename="Atomix",
    modules=[Atomix],
    format=Documenter.HTML(;
        # Only create web pretty-URLs on the CI
        prettyurls=get(ENV, "CI", nothing) == "true",
    ),
    warnonly=:missing_docs,
)

deploydocs(
    repo="github.com/JuliaConcurrent/Atomix.jl",
    devbranch="main",
    push_preview=true,
)

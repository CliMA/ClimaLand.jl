using Documenter, TensorCore

makedocs(;
    modules=[TensorCore],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/JuliaMath/TensorCore.jl/blob/{commit}{path}#L{line}",
    sitename="TensorCore.jl",
    authors="Tim Holy <tim.holy@gmail.com>",
)

deploydocs(;
    repo="github.com/JuliaMath/TensorCore.jl",
    push_preview=true
)

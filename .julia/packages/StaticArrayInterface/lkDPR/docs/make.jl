using StaticArrayInterface
using Documenter

makedocs(;
    modules=[StaticArrayInterface],
    sitename="StaticArrayInterface.jl",
    pages=[
        "StaticArrayInterface.jl: Static Compile-Time Enforced Array Interface Functionality" => "index.md",
        "API" => "api.md"
    ]
)

deploydocs(;
    repo="github.com/JuliaArrays/StaticArrayInterface.jl"
)

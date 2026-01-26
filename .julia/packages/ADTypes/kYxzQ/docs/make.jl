using Documenter
using ADTypes

makedocs(
    sitename = "ADTypes",
    format = Documenter.HTML(),
    modules = [ADTypes]
)

deploydocs(
    repo = "github.com/SciML/ADTypes.jl.git"
)

using Documenter
using NVTX

const ci = get(ENV, "CI", "") == "true"

makedocs(
    sitename = "NVTX.jl",
    format = Documenter.HTML(),
    modules = [NVTX],
    pages = [
        "index.md",
        "api.md",
        "tips.md",
    ]
)

if ci
    deploydocs(
        repo = "github.com/JuliaGPU/NVTX.jl",
        target = "build",
        push_preview = true,
        devbranch = "main",
        forcepush = true,
    )
end

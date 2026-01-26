using Documenter
using ArrayLayouts

DocMeta.setdocmeta!(ArrayLayouts, :DocTestSetup, :(using ArrayLayouts); recursive=true)

makedocs(
    sitename = "ArrayLayouts",
    modules = [ArrayLayouts],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/JuliaLinearAlgebra/ArrayLayouts.jl.git",
)

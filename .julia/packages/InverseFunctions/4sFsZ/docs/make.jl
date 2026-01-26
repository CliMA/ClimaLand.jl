# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [nonstrict] [fixdoctests]
#
# for local builds.

using Documenter
using InverseFunctions

# Doctest setup
DocMeta.setdocmeta!(
    InverseFunctions,
    :DocTestSetup,
    :(using InverseFunctions);
    recursive=true,
)

makedocs(
    sitename = "InverseFunctions",
    modules = [InverseFunctions],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://JuliaMath.github.io/InverseFunctions.jl/stable/"
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    linkcheck = !("nonstrict" in ARGS),
    warnonly = ("nonstrict" in ARGS),
)

deploydocs(
    repo = "github.com/JuliaMath/InverseFunctions.jl.git",
    forcepush = true,
    push_preview = true,
)

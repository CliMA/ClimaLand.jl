using Documenter, DocumenterCitations
using ClimaTimeSteppers
using Literate

tutorial_basedir = "tutorials"
tutorial_basedir_from_here = joinpath(@__DIR__, "src", tutorial_basedir)

jl_files_in_basedir = filter(endswith(".jl"), readdir(tutorial_basedir_from_here))

println("Building literate tutorials...")
generated_tutorials = String[]
for filename in jl_files_in_basedir
    Literate.markdown(
        joinpath(tutorial_basedir_from_here, filename),
        tutorial_basedir_from_here;
        execute = true,
        flavor = Literate.CommonMarkFlavor(),
    )
    push!(generated_tutorials, joinpath(tutorial_basedir, replace(filename, ".jl" => ".md")))
end

# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "nul"

bib = CitationBibliography(joinpath(@__DIR__, "refs.bib"))

#! format: off
pages = [
    "index.md",
    "Algorithm Formulations" => [
        "ODE Solvers" => "algorithm_formulations/ode_solvers.md",
        "Newtons Method" => "algorithm_formulations/newtons_method.md",
        "Rosenbrock Method" => "algorithm_formulations/rosenbrock.md",
        "Old LSRK Formulations" => "algorithm_formulations/lsrk.md",
        "Old MRRK Formulations" => "algorithm_formulations/mrrk.md",
    ],
    "Test problems" => [
        "test_problems/index.md",
    ],
    "Tutorials" => generated_tutorials,
    "API docs" => [
        "ODE Solvers" => "api/ode_solvers.md",
        "Newtons Method" => "api/newtons_method.md",
        "Callbacks" => "api/callbacks.md",
    ],
    # "Algorithm comparisons" => "algo_comparisons.md", # TODO: fill out
    "Developer docs" => [
        "Types" => "dev/types.md",
        "Report generator" => "dev/report_gen.md",
    ],
    "references.md",
]
#! format: on

mathengine = MathJax(Dict(:TeX => Dict(:equationNumbers => Dict(:autoNumber => "AMS"), :Macros => Dict())))

format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true", mathengine = mathengine, collapselevel = 1)

makedocs(;
    plugins = [bib],
    sitename = "ClimaTimeSteppers",
    format = format,
    modules = [ClimaTimeSteppers],
    checkdocs = :exports,
    clean = true,
    pages = pages,
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/CliMA/ClimaTimeSteppers.jl.git",
    target = "build",
    push_preview = true,
    devbranch = "main",
    forcepush = true,
)

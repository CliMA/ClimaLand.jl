# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "nul"
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using Distributed
@everywhere using Documenter
@everywhere using Literate
@everywhere using ClimaLand
include("pages_helper.jl")
include("list_tutorials.jl")

@everywhere const clima_dir = dirname(dirname(pathof(ClimaLand)));
@everywhere source_dir = joinpath(@__DIR__, "src")
@everywhere GENERATED_DIR = joinpath(source_dir, "generated") # generated files directory
rm(GENERATED_DIR, force = true, recursive = true)
mkpath(GENERATED_DIR)
@everywhere function generate_tutorial(tutorials_dir, tutorial)
    rpath = relpath(dirname(tutorial), tutorials_dir)
    rpath = rpath == "." ? "" : rpath
    gen_dir = joinpath(GENERATED_DIR, rpath)
    mkpath(gen_dir)

    cd(gen_dir) do
        # change the Edit on GitHub link:
        path = relpath(clima_dir, pwd())
        content = """
                # ```@meta
                    # EditURL = "https://github.com/CliMA/ClimaLand.jl/$(path)"
                    # ```
                """
        mdpre(str) = content * str
        input = abspath(tutorial)
        Literate.markdown(
            input;
            execute = true,
            documenter = false,
            preprocess = mdpre,
        )
    end
end
tutorials_jl = flatten_to_array_of_strings(get_second(tutorials))
println("Building literate tutorials...")
tutorials_dir = joinpath(@__DIR__, "src", "tutorials")
tutorials_jl = map(x -> joinpath(tutorials_dir, x), tutorials_jl)
pmap(t -> generate_tutorial(tutorials_dir, t), tutorials_jl)

# update list of rendered markdown tutorial output for mkdocs
ext_jl2md(x) = joinpath(basename(GENERATED_DIR), replace(x, ".jl" => ".md"))
tutorials = transform_second(x -> ext_jl2md(x), tutorials)
include("list_of_apis.jl")
include("list_standalone_models.jl")
include("list_diagnostics.jl")
pages = Any[
    "Home" => "index.md",
    "Running your first simulation" => "getting_started.md",
    "Tutorials" => tutorials,
    "Model Equations" => standalone_models,
    "Additional resources" => [
        "Model structure" => "model_structure.md",
        "Repository structure" => "repo_structure.md",
        "Analyzing model output" => [
            "Diagnostics" => diagnostics,
            "Leaderboard" => "leaderboard/leaderboard.md",
        ],
        "Running on GPU or with MPI" => "architectures.md",
        "Calibration" => "calibration.md",
        "Restarting a simulation" => "restarts.md",
        "Software utilities" => [
            "ITime type" => "itime.md",
            "Shared utilities" => "shared_utilities.md",
        ],
        "Physical units" => "physical_units.md",
        "Julia background" => "julia.md",
    ],
    "APIs" => apis,
    "Contributor guide" => "contributing.md",
]

mathengine = MathJax(
    Dict(
        :TeX => Dict(
            :equationNumbers => Dict(:autoNumber => "AMS"),
            :Macros => Dict(),
        ),
    ),
)

format = Documenter.HTML(
    prettyurls = !isempty(get(ENV, "CI", "")),
    collapselevel = 1,
    mathengine = mathengine,
    ansicolor = true,
    size_threshold = 500_000,  # 500 KiB instead of default 200 KiB
    size_threshold_warn = 300_000,  # 300 KiB instead of default 100 KiB
)

makedocs(
    sitename = "ClimaLand.jl",
    authors = "Clima Land Model Team",
    format = format,
    pages = pages,
    checkdocs = :exports,
    doctest = true,
    warnonly = [:missing_docs, :footnote, :cross_references],
    clean = true,
    modules = [ClimaLand],
)

deploydocs(
    repo = "github.com/CliMA/ClimaLand.jl.git",
    target = "build",
    push_preview = true,
    devbranch = "main",
    forcepush = true,
)

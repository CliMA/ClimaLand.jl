# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "nul"
push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using Distributed
@everywhere using Documenter
@everywhere using Literate
using ClimaLand
include("pages_helper.jl")
tutorials = [
    "For model developers" => [
        "Intro to standalone models" => "standalone/Usage/model_tutorial.jl",
        "Intro to multi-component models" => "standalone/Usage/LSM_single_column_tutorial.jl",
        "Intro to ClimaLand Domains" => "standalone/Usage/domain_tutorial.jl",
        "Intro to forced site-level runs" => "shared_utilities/driver_tutorial.jl",
    ],
    "Running simulations" => [
        "Bucket LSM" => [
            "standalone/Bucket/bucket_tutorial.jl",
            "standalone/Bucket/coupled_bucket.jl",
        ],
        "Soil modeling" => [
            "Boundary conditions" => "standalone/Soil/boundary_conditions.jl",
            "Richards Equation" => "standalone/Soil/richards_equation.jl",
            "Energy and Hydrology" => "standalone/Soil/soil_energy_hydrology.jl",
            "Phase Changes" => "standalone/Soil/freezing_front.jl",
            "Layered Soil" => "standalone/Soil/layered_soil.jl",
            "Coarse Sand Evaporation" => "standalone/Soil/evaporation.jl",
            "Gilat Loess Evaporation" => "standalone/Soil/evaporation_gilat_loess.jl",
            "Bare soil site" => "standalone/Soil/sublimation.jl",
        ],
        "Canopy modeling" => [
            "Standalone Canopy" => "standalone/Canopy/canopy_tutorial.jl",
        ],
        "Integrated soil+canopy modeling" => [
            "Coupled Canopy and Soil" => "integrated/soil_canopy_tutorial.jl",
        ],
        "Bucket LSM" => [
            "standalone/Bucket/bucket_tutorial.jl",
            "standalone/Bucket/coupled_bucket.jl",
        ],
        "Snow Modeling" => [
            "standalone/Snow/base_tutorial.jl",
            "standalone/Snow/data_tutorial.jl",
        ],
    ],
]
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
tutorials_dir = joinpath(@__DIR__, "tutorials")
tutorials_jl = map(x -> joinpath(tutorials_dir, x), tutorials_jl)
pmap(t -> generate_tutorial(tutorials_dir, t), tutorials_jl)

# update list of rendered markdown tutorial output for mkdocs
ext_jl2md(x) = joinpath(basename(GENERATED_DIR), replace(x, ".jl" => ".md"))
tutorials = transform_second(x -> ext_jl2md(x), tutorials)
include("list_of_apis.jl")
include("list_standalone_models.jl")
pages = Any[
    "Home" => "index.md",
    "APIs" => apis,
    "Contribution guide" => "Contributing.md",
    "Tutorials" => tutorials,
    "Repository structure" => "folderstructure.md",
    "Standalone models" => standalone_models,
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
)

makedocs(
    sitename = "ClimaLand.jl",
    authors = "Clima Land Model Team",
    format = format,
    pages = pages,
    checkdocs = :exports,
    doctest = true,
    warnonly = true,
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

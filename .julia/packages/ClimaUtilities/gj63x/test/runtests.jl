using SafeTestsets
using Pkg

# Download test-only artifacts
#
# (Currently not natively supported by Julia)
artifacts_toml = joinpath(@__DIR__, "Artifacts.toml")
artifacts = Pkg.Artifacts.select_downloadable_artifacts(artifacts_toml)
for name in keys(artifacts)
    Pkg.Artifacts.ensure_artifact_installed(
        name,
        artifacts[name],
        artifacts_toml,
    )
end

# Performance and code quality tests
@safetestset "Quality Assurance tests" begin
    include("quality_assurance.jl")
end

# Unit tests
@safetestset "Utils tests" begin
    include("utils.jl")
end

@safetestset "OutputPathGenerator tests" begin
    include("output_path_generator.jl")
end

@safetestset "TimeManager tests" begin
    include("timemanager.jl")
    include("itime.jl")
end

@safetestset "DataStructures tests" begin
    include("data_structures.jl")
end

@safetestset "FileReaders tests" begin
    include("file_readers.jl")
end

@safetestset "Regridders tests" begin
    include("regridders.jl")
end

@safetestset "DataHandling tests" begin
    include("data_handling.jl")
end

@safetestset "SpaceVaryingInputs tests" begin
    include("space_varying_inputs.jl")
end

@safetestset "TimeVaryingInputs tests" begin
    include("time_varying_inputs.jl")
end

@safetestset "TimeVaryingInputs23D tests" begin
    include("time_varying_inputs23.jl")
end

@safetestset "TimeVaryingInputs tests LinearPeriodFilling" begin
    include("time_varying_inputs_linearperiodfilling.jl")
end

@safetestset "ClimaArtifacts tests" begin
    include("clima_artifacts.jl")
end

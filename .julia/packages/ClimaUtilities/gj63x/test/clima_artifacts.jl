using Artifacts
using LazyArtifacts
using Test

import ClimaUtilities.ClimaArtifacts

# Test with a non lazy-artifact (and without ClimaComms)
expected_path2 = artifact"laskar2004"

@testset "Non-lazy artifact, context: (without ClimaComms)" begin
    @test ClimaArtifacts.@clima_artifact("laskar2004") == expected_path2

    # Test with name as a variable
    artifact_name = "laskar2004"

    @test ClimaArtifacts.@clima_artifact(artifact_name) == expected_path2
end

import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends

const context = ClimaComms.context()
ClimaComms.init(context)

let_filesystem_catch_up() = context isa ClimaComms.MPICommsContext && sleep(0.2)

expected_path = artifact"socrates"

# Remove the artifact, so that we test that we are downloading it
ClimaComms.iamroot(context) &&
    Base.Filesystem.rm(expected_path, force = true, recursive = true)
ClimaComms.barrier(context)
let_filesystem_catch_up()
@info "Removed artifact"

@testset "Lazy artifact with context" begin
    @test_throws ErrorException @macroexpand ClimaArtifacts.@clima_artifact(
        "socrates"
    )

    @test ClimaArtifacts.@clima_artifact("socrates", context) == expected_path

    # Test with name as a variable
    ClimaComms.iamroot(context) &&
        Base.Filesystem.rm(expected_path, force = true, recursive = true)
    ClimaComms.barrier(context)
    let_filesystem_catch_up()
    artifact_name = "socrates"

    @test ClimaArtifacts.@clima_artifact(artifact_name, context) ==
          expected_path

    @test_throws ErrorException ClimaArtifacts.@clima_artifact(artifact_name)
end

@testset "Accessed artifacts" begin
    @test ClimaArtifacts.accessed_artifacts() == Set(["socrates", "laskar2004"])
end

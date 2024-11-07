"""
This module provides functions for accessing test data files used solely in the
ClimaLand test suite. Each function provides the path to a specific test data
file which may then be read in and used in the relevant test scripts.
"""

import Pkg, Artifacts
import ClimaUtilities.ClimaArtifacts: @clima_artifact

# Download test-only artifacts
#
# (Currently not natively supported by Julia)
artifacts_toml = joinpath(@__DIR__, "Artifacts.toml")
artifacts = Artifacts.select_downloadable_artifacts(artifacts_toml)
for name in keys(artifacts)
    Pkg.Artifacts.ensure_artifact_installed(
        name,
        artifacts[name],
        artifacts_toml,
    )
end

"""
    twostr_test_data_path(; context = nothing)

Returns the file path for data used to test the ClimaLand TwoStream model
implementation. This test data is produced from the PySellersTwoStream
implementation by T. Quaife, and the test data is used to check that the
outputs of the two models are equal.
Quaife, T. 2016: PySellersTwoStream, available at:
https://github.com/tquaife/pySellersTwoStream
"""
function twostr_test_data_path(; context = nothing)
    return @clima_artifact("twostr_test", context)
end

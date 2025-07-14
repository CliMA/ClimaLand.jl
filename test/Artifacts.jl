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

"""
    water_conservation_test_data_path(; context = nothing)

Returns the filepaths for data from two simulations of ClimaLand.Soil.RichardsModel;
these were carried out with a very small timestep with an explicit timestepper
and are used as ground truth for solutions using an implicit timestepper.

Experiment details are in `experiments/standalone/Soil/water_conservation.jl`.
"""
function water_conservation_test_data_path(; context = nothing)
    folder_path = @clima_artifact("water_conservation_test", context)
    flux_datapath = joinpath(folder_path, "ref_soln_flux.csv")
    dirichlet_datapath = joinpath(folder_path, "ref_soln_dirichlet.csv")
    return flux_datapath, dirichlet_datapath
end


"""
    pmodel_unittests_path(; context = nothing)

Returns the filepaths for input-output pairs for the P-model. This data is generated from the R
implementation of the P-model and is used to test the ClimaLand P-model.
    
Experiment details are in `test/standalone/Vegetation/test_pmodel.jl`.
"""
function pmodel_unittests_path(; context = nothing)
    return @clima_artifact("pmodel_unittests", context)
end

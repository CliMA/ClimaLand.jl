using ClimaLand
using ClimaLand.Artifacts
using DelimitedFiles
using Logging
using Base.Threads

using CUDA
using ClimaComms
ClimaComms.@import_required_backends()

script = joinpath(pkgdir(ClimaLand), "experiments/calibration/test_calibration.jl")
include(script)  # to load dependencies

# Get all site IDs from the metadata file
fluxnet2015_data_path = ClimaLand.Artifacts.fluxnet2015_data_path()
fluxnet2015_metadata_path =
        joinpath(fluxnet2015_data_path, "metadata_DD_clean.csv")

data = readdlm(fluxnet2015_metadata_path, ',', header=true)

header = data[2]                  # column names
table = data[1]                   # the actual rows

site_id_col = findfirst(==("site_id"), header)
site_ids = table[:, site_id_col]  # vector of all site IDs
unique_site_ids = unique(site_ids)

# Run the script for each unique site ID
Threads.@threads for i in 1:length(unique_site_ids) 
    site_ID = unique_site_ids[i]

    @info "Running $script with file: $site_ID."
    try
        # run(`julia --project=.buildkite $script $site_ID`)
        calibrate_at_site(site_ID)
        @info "Success on $site_ID"
    catch e
        @info "Failed on $site_ID: $e"
        continue
    end
end
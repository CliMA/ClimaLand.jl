using ClimaLand
using ClimaLand.Artifacts
using DelimitedFiles
using Logging

script = joinpath(pkgdir(ClimaLand), "experiments/calibration/test_calibration.jl")

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
for site_ID in unique_site_ids
    println("Running $script with file: $site_ID")
    try
        run(`julia --project=.buildkite $script $site_ID`)
        println("Success on $site_ID")
    catch e
        println("Failed on $site_ID: $e")
        continue
    end
end
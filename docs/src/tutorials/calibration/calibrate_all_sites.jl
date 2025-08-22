using ClimaLand
using ClimaLand.Artifacts
using DelimitedFiles

# folder = "path/to/folder"
script = "test_calibration.jl"

fluxnet2015_data_path = ClimaLand.Artifacts.fluxnet2015_data_path()
fluxnet2015_metadata_path =
        joinpath(fluxnet2015_data_path, "metadata_DD_clean.csv")


data = readdlm(fluxnet2015_metadata_path, ',', header=true)

header = data[2]                  # column names
table = data[1]                   # the actual rows

site_id_col = findfirst(==("SITE_ID"), header)
site_ids = table[:, site_id_col]  # vector of all site IDs
unique_site_ids = unique(site_ids)

for f in readdir(folder; join=true)
    try
        println("Running $script with file: $f")
        run(`julia $script $f`)
        println("✅ Success on $f")
    catch e
        println("❌ Failed on $f: $e")
        continue
    end
end
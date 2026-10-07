# Creates `clm_crop_fraction.nc`: the fraction of the land in the crop land unit in the
# CLM5 surface data for the year 2000 (0.9°×1.25°, `PCT_CROP`, from LUH2), NaN over
# ocean. A temporary copy until this becomes a ClimaArtifacts artifact.
#
#     julia --project=.buildkite artifacts/clm_crop_fraction/create_clm_crop_fraction.jl [surfdata]
#
# `surfdata` is a local copy of the surface data file, downloaded if not given.

import Downloads
import NCDatasets
import SHA

const URL = "https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/lnd/clm2/surfdata_map/surfdata_0.9x1.25_16pfts__CMIP6_simyr2000_c170616.nc"
const SHA256 = "b1062bac7a907e5841ec4217cc29a9159758c40c385a4248125a958caa1891c8"

source = isempty(ARGS) ? Downloads.download(URL) : only(ARGS)
@assert bytes2hex(open(SHA.sha256, source)) == SHA256

ds = NCDatasets.NCDataset(source)
lon = Array(ds["LONGXY"][:, 1])
lat = Array(ds["LATIXY"][1, :])
pct_crop = Array(ds["PCT_CROP"][:, :])
land_fraction = Array(ds["LANDFRAC_PFT"][:, :])
close(ds)
crop_fraction = @. ifelse(land_fraction > 0, pct_crop / 100, NaN)

NCDatasets.NCDataset(joinpath(@__DIR__, "clm_crop_fraction.nc"), "c") do out
    out.attrib["title"] = "Crop fraction of the land, CLM5 surface data, year 2000"
    out.attrib["source"] = URL
    NCDatasets.defVar(
        out,
        "lon",
        lon,
        ("lon",);
        attrib = ["units" => "degrees_east"],
    )
    NCDatasets.defVar(
        out,
        "lat",
        lat,
        ("lat",);
        attrib = ["units" => "degrees_north"],
    )
    NCDatasets.defVar(
        out,
        "crop_fraction",
        Float32.(crop_fraction),
        ("lon", "lat");
        deflatelevel = 9,
        attrib = [
            "units" => "fraction",
            "long_name" => "Fraction of the land in the crop land unit",
        ],
    )
end

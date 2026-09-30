# Creates `xusaatchi_woody_carbon_2deg.nc`: the 2000-2019 mean live woody carbon of
# Xu et al. (2021), https://doi.org/10.5281/zenodo.4161694, as distributed by ILAMB at
# 0.5°, averaged in 2° blocks and converted to kg C m^-2. A temporary copy until this
# becomes a ClimaArtifacts artifact.
#
#     julia --project=.buildkite artifacts/prognostic_carbon/create_xusaatchi_woody_carbon.jl

import Downloads
import NCDatasets
import SHA
import Statistics: mean

const URL = "https://www.ilamb.org/ILAMB-Data/DATA/biomass/XuSaatchi2021/XuSaatchi.nc"
const SHA256 = "ab93cb57ee2a9e35fb0fff2e64da2f03a6a83358ef9936f52e51f929f143cc36"
const BLOCK = 4 # 0.5° cells per 2° block

source = Downloads.download(URL)
@assert bytes2hex(open(SHA.sha256, source)) == SHA256

ds = NCDatasets.NCDataset(source)
lon = Array(ds["lon"][:])
lat = Array(ds["lat"][:])
# Mg C ha^-1 to kg C m^-2, averaged over the years
biomass = 0.1 .* ds["biomass"].var[:, :, :]
close(ds)
valid(x) = isfinite(x) && x >= 0
time_mean = map(CartesianIndices(size(biomass)[1:2])) do I
    x = filter(valid, biomass[I, :])
    isempty(x) ? NaN : mean(x)
end
# Latitudes ascending
order = sortperm(lat)
lat, time_mean = lat[order], time_mean[:, order]

block_mean(x) = (v = filter(isfinite, x); isempty(v) ? NaN : mean(v))
nlon, nlat = length(lon) ÷ BLOCK, length(lat) ÷ BLOCK
blocks(i) = ((i - 1) * BLOCK + 1):(i * BLOCK)
woody_carbon =
    [block_mean(time_mean[blocks(i), blocks(j)]) for i in 1:nlon, j in 1:nlat]
block_lon = [mean(lon[blocks(i)]) for i in 1:nlon]
block_lat = [mean(lat[blocks(j)]) for j in 1:nlat]

NCDatasets.NCDataset(
    joinpath(@__DIR__, "xusaatchi_woody_carbon_2deg.nc"),
    "c",
) do out
    out.attrib["title"] = "Live woody carbon, 2000-2019 mean, 2 degree blocks"
    out.attrib["source"] = URL
    out.attrib["references"] = "Xu et al. (2021), https://doi.org/10.5281/zenodo.4161694"
    NCDatasets.defVar(
        out,
        "lon",
        block_lon,
        ("lon",);
        attrib = ["units" => "degrees_east"],
    )
    NCDatasets.defVar(
        out,
        "lat",
        block_lat,
        ("lat",);
        attrib = ["units" => "degrees_north"],
    )
    NCDatasets.defVar(
        out,
        "woody_carbon",
        Float32.(woody_carbon),
        ("lon", "lat");
        deflatelevel = 9,
        attrib = [
            "units" => "kg m-2",
            "long_name" => "Carbon in live woody vegetation",
        ],
    )
end

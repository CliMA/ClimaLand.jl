# Creates `clm_tree_share.nc`: the tree share of the natural vegetation in the CLM5
# surface data for the year 2000 (0.9°×1.25°), the percent cover of the tree PFTs
# (1-8) over that of all vegetated PFTs (1-14, excluding bare ground), NaN where there
# is no vegetation. A temporary copy until this becomes a ClimaArtifacts artifact.
#
#     julia --project=.buildkite artifacts/clm_tree_share/create_clm_tree_share.jl [surfdata]
#
# `surfdata` is a local copy of the surface data file, downloaded if not given.

import Downloads
import NCDatasets
import SHA

const URL = "https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/lnd/clm2/surfdata_map/surfdata_0.9x1.25_16pfts__CMIP6_simyr2000_c170616.nc"
const SHA256 = "b1062bac7a907e5841ec4217cc29a9159758c40c385a4248125a958caa1891c8"
# Indices into the natpft dimension (0-14) of PCT_NAT_PFT
const TREE_PFTS = 2:9
const VEGETATED_PFTS = 2:15

source = isempty(ARGS) ? Downloads.download(URL) : only(ARGS)
@assert bytes2hex(open(SHA.sha256, source)) == SHA256

ds = NCDatasets.NCDataset(source)
lon = Array(ds["LONGXY"][:, 1])
lat = Array(ds["LATIXY"][1, :])
pct = Array(ds["PCT_NAT_PFT"][:, :, :])
close(ds)
trees = dropdims(sum(pct[:, :, TREE_PFTS]; dims = 3); dims = 3)
vegetated = dropdims(sum(pct[:, :, VEGETATED_PFTS]; dims = 3); dims = 3)
tree_share = @. ifelse(vegetated > 0, trees / vegetated, NaN)

NCDatasets.NCDataset(joinpath(@__DIR__, "clm_tree_share.nc"), "c") do out
    out.attrib["title"] = "Tree share of the natural vegetation, CLM5 surface data, year 2000"
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
        "tree_share",
        Float32.(tree_share),
        ("lon", "lat");
        deflatelevel = 9,
        attrib = [
            "units" => "fraction",
            "long_name" => "Tree share of the vegetated natural cover",
        ],
    )
end

using NCDatasets
using ClimaArtifactsHelper
using Statistics

if length(ARGS) < 1
    @error("Please provide the local path to ClimaLand simulation diagnostics.")
else
    filedir = ARGS[1]
end
outdir = "rosetta_spunup_ic"
if isdir(outdir)
    @warn "$outdir already exists. Content will end up in the artifact and may be overwritten."
    @warn "Abort this calculation, unless you know what you are doing."
else
    mkdir(outdir)
end

# First, we read in the output from a ClimaLand long run with the snowy land model.
# This was forced with ERA5 data for 20 years using the ROSETTA soil parameters.

attrib_swc = (;
                   vartitle = "Volumetric fraction of water",
                   varunits = "m^3/m^3",
                   varname = "swc",
                   )
attrib_si = (;
                   vartitle = "Volumetric fraction of ice",
                   varunits = "m^3/m^3",
                   varname = "si",
                   )
attrib_sie = (;
                   vartitle = "Soil volumetric internal energy",
                   varunits = "J/m^3",
                   varname = "sie",
              )
attrib_swe = (;
                   vartitle = "Snow water equivalent",
                   varunits = "m",
                   varname = "swe",
              )
	      attrib_lwp = (;
                   vartitle = "Leaf water potential",
                   varunits = "m",
                   varname = "lwp",
              )

var_3d_attribs = [attrib_swc, attrib_si, attrib_sie]
var_2d_attribs = [attrib_swe, attrib_lwp]
function replace_nan_with_mean!(x)
    nan_mask = isnan.(x)
    nonnan_mean = mean(x[.~nan_mask])
    x[nan_mask] .= nonnan_mean
    return nothing
end

outfilepath = joinpath(outdir, "rosetta_spunup_ic_20years.nc")
ds = NCDataset(outfilepath, "c")    

for i in 1:length(var_3d_attribs)
    var_attrib = var_3d_attribs[i]
    (vartitle, varunits, varname) = var_attrib
    data = NCDataset(joinpath(filedir, "$(varname)_1M_average.nc"))
    z = data["z"][:]
    lat = data["lat"][:]
    lon = data["lon"][:]
        if i == 1
            # Define and set values for dimensions (lon, lat, z)
            defDim(ds, "lon", length(lon))
            defDim(ds, "lat", length(lat))
            defDim(ds, "z", length(z))
            la = defVar(ds, "lat", Float32, ("lat",))
            lo = defVar(ds, "lon", Float32, ("lon",))
            zv = defVar(ds, "z", Float32, ("z",))
            la.attrib["units"] = "degrees_north"
            la.attrib["standard_name"] = "latitude"
            lo.attrib["standard_name"] = "longitude"
            lo.attrib["units"] = "degrees_east"
            zv.attrib["standard_name"] = "depth"
            zv.attrib["units"] = "m"
            la[:] = lat
            lo[:] = lon
            zv[:] = z
        end
    var = defVar(ds, varname, Float32, ("lon", "lat", "z"))
    var.attrib["units"] = varunits
    var.attrib["longname"]= vartitle
    var.attrib["varname"] = varname
    field = data[varname][end,:,:,:]; # last time element
    close(data)
    replace_nan_with_mean!(field)
    var[:, :, :] = field[:,:,:]
end
for i in 1:length(var_2d_attribs)
    var_attrib = var_2d_attribs[i]
    (vartitle, varunits, varname) = var_attrib
    data = NCDataset(joinpath(filedir, "$(varname)_1M_average.nc"))
    lat = data["lat"][:]
    lon = data["lon"][:]
    var = defVar(ds, varname, Float32, ("lon", "lat"))
    var.attrib["units"] = varunits
    var.attrib["longname"]= vartitle
    var.attrib["varname"] = varname
    field = data[varname][end,:,:]; # last time element
    close(data)
    replace_nan_with_mean!(field)
    var[:, :] = field[:,:]
end

close(ds)

#create_artifact_guided(outdir; artifact_name = basename(@__DIR__))
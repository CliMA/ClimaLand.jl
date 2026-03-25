using NCDatasets
using ClimaLand
using Dates
start_date = DateTime(1980)
stop_date = DateTime(2020)
era5_ncdata_path = ClimaLand.Artifacts.find_era5_year_paths(
    start_date,
    stop_date;
)
field = zeros((360,181))
count = 0
for file in era5_ncdata_path
    @show file
    data = NCDataset(file)
    if file == era5_ncdata_path[1]
        lat = data["lat"][:]
        lon = data["lon"][:]
    end
    field .+= sum(data["t2m"][:,:,:], dims = 3)[:,:,1]/length(data["time"][:])
    count +=1
    close(data)
end
field = field ./ count

outfilepath = joinpath("mean_sfc_air_temp.nc")
ds = NCDataset(outfilepath, "c")
defDim(ds, "lon", length(lon))
defDim(ds, "lat", length(lat))
la = defVar(ds, "lat", Float32, ("lat",))
lo = defVar(ds, "lon", Float32, ("lon",))
la.attrib["units"] = "degrees_north"
la.attrib["standard_name"] = "latitude"
lo.attrib["standard_name"] = "longitude"
lo.attrib["units"] = "degrees_east"
la[:] = lat
lo[:] = lon
attrib_lwp = (;
              vartitle = "Mean air temp",
              varunits = "K",
              varname = "tair",
              )
var = defVar(ds, varname, Float32, ("lon", "lat"))
var.attrib["units"] = varunits
var.attrib["longname"]= vartitle
var.attrib["varname"] = varname
var[:, :] = field[:,:]
close(ds)

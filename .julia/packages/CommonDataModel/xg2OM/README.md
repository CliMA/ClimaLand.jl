[![Build Status](https://github.com/JuliaGeo/CommonDataModel.jl/workflows/CI/badge.svg)](https://github.com/JuliaGeo/CommonDataModel.jl/actions)
[![codecov](https://codecov.io/github/JuliaGeo/CommonDataModel.jl/graph/badge.svg?token=TNU4HSPelE)](https://codecov.io/github/JuliaGeo/CommonDataModel.jl)
[![documentation stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliageo.github.io/CommonDataModel.jl/stable/)
[![documentation dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliageo.github.io/CommonDataModel.jl/dev/)


This package contains abstracts type definition for loading and manipulating GRIB, NetCDF, geoTiff and Zarr files. This package aims to follow the [Common Data Model](https://docs.unidata.ucar.edu/netcdf-c/current/netcdf_data_model.html) and the [CF (climate and forecast models) Metadata Conventions](https://cfconventions.org/).


| Format  |      Package | read support | write support |
|---------|--------------|:------------:|:-------------:|
| NetCDF  | [`NCDatasets`](https://github.com/Alexander-Barth/NCDatasets.jl)     |            ✔ |             ✔ |
| OPeNDAP | [`NCDatasets`](https://github.com/Alexander-Barth/NCDatasets.jl)     |            ✔ |             - |
| GRIB    | [`GRIBDatasets`](https://github.com/JuliaGeo/GRIBDatasets.jl)        |            ✔ |             - |
| geoTIFF | [`TIFFDatasets`](https://github.com/Alexander-Barth/TIFFDatasets.jl) |            ✔ |             - |
| Zarr    | [`ZarrDatasets`](https://github.com/JuliaGeo/ZarrDatasets.jl)        |            ✔ |             ✔ |
|         | [`MetopDatasets`](https://github.com/eumetsat/MetopDatasets.jl)        |            ✔ |             - |


Features include:
* query and edit metadata of arrays and datasets
* virtually concatenating multiple files along a given dimension and merging virtually different datasets
* create a virtual subset (`view`) by indices or by values of coordinate variables (`CommonDataModel.select`, `CommonDataModel.@select`)
* group, map and reduce a variable (`CommonDataModel.groupby`, `CommonDataModel.@groupby`) and rolling reductions like running means `CommonDataModel.rolling`)




Here is minimal example for loading files using `CommonDataModel`:

``` julia
import CommonDataModel as CDM
import SomeDatasets # where SomeDatasets is either GRIBDatasets, NCDatasets, ZarrDatasets,...

ds = SomeDatasets.Dataset("file_name")

# ntime is the number of time instances
ntime = ds.dim["time"] # or CDM.dims(ds)["time"]

# create an array-like structure v corresponding to variable temperature
v = ds["temperature"]

# load a subset
subdata = v[10:30,30:5:end]

# load all data
data = v[:,:]

# load a global attribute
title = ds.attrib["title"]  # or CDM.attribs(ds)["title"]
close(ds)
```

Most users would typically import [`GRIBDatasets`](https://github.com/JuliaGeo/GRIBDatasets.jl), [`NCDatasets`](https://github.com/Alexander-Barth/NCDatasets.jl)... directly and not `CommonDataModel`.

# File conversions

By implementing a common interface, files can be converted from one format to another using the `write` function.
For example GRIB files can be converted to NetCDF (or Zarr) files:

```julia
using NCDatasets # or ZarrDatasets
using GRIBDatasets
using Downloads: download

grib_file = download("https://github.com/JuliaGeo/GRIBDatasets.jl/raw/98356af026ea39a5ec0b5e64e4289105492321f8/test/sample-data/era5-levels-members.grib")
netcdf_file = "test.nc"
NCDataset(netcdf_file,"c") do ds
   write(ds,GRIBDataset(grib_file))
end
```

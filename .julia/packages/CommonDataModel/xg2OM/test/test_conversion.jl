import NCDatasets
import GRIBDatasets
import CommonDataModel as CDM
using Test
using DataStructures

datadir = joinpath(dirname(pathof(GRIBDatasets)),"..","test","sample-data")
filename = joinpath(datadir,"era5-levels-members.grib")
ds = GRIBDatasets.Dataset(filename)

io = IOBuffer()
show(io,ds)
out = String(take!(io))
@test occursin("Global attributes",out)
@test occursin("CF-",out)

# dimension

@test CDM.dims(ds)["lon"] == 120
@test CDM.dim(ds,"lon") == 120
@test "lon" in CDM.dimnames(ds)


tmp_filename = tempname()
NCDatasets.write(tmp_filename,ds)

@test isfile(tmp_filename)

dsnc = NCDatasets.Dataset(tmp_filename)
@test ds["number"][:] == dsnc["number"][:]

@test OrderedDict(CDM.dims(ds)) == OrderedDict(CDM.dims(dsnc))
@test OrderedDict(CDM.attribs(ds)) == OrderedDict(CDM.attribs(dsnc))
@test OrderedDict(CDM.groups(ds)) == OrderedDict(CDM.groups(dsnc))


close(dsnc)


# copy individual variables
tmp_filename = tempname()
dsnc = NCDatasets.Dataset(tmp_filename,"c")
gribv = ds["lon"]
ncv = CDM.defVar(dsnc,gribv)
@test CDM.name(ncv) == CDM.name(gribv)
@test ncv[:] == gribv[:]
@test CDM.attrib(ncv,"units") == CDM.attrib(gribv,"units")
@test collect(CDM.dimnames(ncv)) == collect(CDM.dimnames(gribv))
close(dsnc)

close(ds)

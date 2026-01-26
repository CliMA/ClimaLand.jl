using Test
using Dates
using Printf
using DataStructures
using CFTime

sz = (4,5)
filename = tempname()

#TDS = NCDatasets.NCDataset
TDS = MemoryDataset

ds = TDS(filename,"c")
data = rand(1:10,sz)
v = defVar(ds,"data",data,("lon","lat"))



filename2 = tempname()
ds2 = TDS(filename2,"c")
ds2["new_data"] = ds["data"]
@test ds2["new_data"][:,:] == ds["data"][:,:]
close(ds2)


filename2 = tempname()
ds2 = TDS(filename2,"c")
ds2["new_data"] = view(ds["data"],1:2,1:2)
@test ds2["new_data"][:,:] == ds["data"][1:2,1:2]
close(ds2)

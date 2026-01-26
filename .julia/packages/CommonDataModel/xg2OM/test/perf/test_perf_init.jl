# create the test file

using NCDatasets
using Dates

sz = (360,180)
time = DateTime(1980,1,1):Day(1):DateTime(2010,1,1);
varname = "data"

fname = expanduser("~/sample_perf2.nc")

isfile(fname) && rm(fname)

NCDataset(fname,"c") do ds
    defVar(ds,"lon",1:sz[1],("lon",))
    defVar(ds,"lat",1:sz[2],("lat",))
    defVar(ds,"time",time,("time",))
    ncv = defVar(ds,"data",Float32,("lon","lat","time"),attrib=Dict("foo" => "bar"))

    for n = 1:length(time)
        ncv[:,:,n] = randn(Float32,sz...) .+ 100
    end
end

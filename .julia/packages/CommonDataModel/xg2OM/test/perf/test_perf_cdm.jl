# Benchmark to be run on Linux as root

using BenchmarkTools
using NCDatasets
using Dates
using CommonDataModel: @groupby
using CommonDataModel

fname = expanduser("~/sample_perf2.nc")
ds = NCDataset(fname)

data_f64 = Float64.(ds[:data][:,:,:])

println("runtime of mean")
gm = @btime begin
    write("/proc/sys/vm/drop_caches","3")
    mean(@groupby(ds[:data],Dates.Month(time)))[:,:,:];
end

println("runtime of std")

# Welford
gs = @btime begin
    write("/proc/sys/vm/drop_caches","3")
    std(@groupby(ds[:data],Dates.Month(time)))[:,:,:];
end

println("accuracy")

mean_ref = cat(
        [mean(data_f64[:,:,findall(Dates.month.(ds[:time][:]) .== m)],dims=3)
         for m in 1:12]...,dims=3);

std_ref = cat(
        [std(data_f64[:,:,findall(Dates.month.(ds[:time][:]) .== m)],dims=3)
         for m in 1:12]...,dims=3);

@show sqrt(mean((gm - mean_ref).^2))
@show sqrt(mean((gs - std_ref).^2))

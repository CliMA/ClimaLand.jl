using Test
using Dates

using CommonDataModel:
    GroupMapping,
    MClass,
    OverlappingGroupMapping,
    dataindices,
    rolling,
    groupindices,
    ndata,
    ngroups

using Statistics


#include("memory_dataset.jl");

function test_mapping(gmap,coord,rolling_classes)
    for k = 1:ndata(gmap)
        for ku in groupindices(gmap,k)
            @test time[k] in rolling_classes[ku]
        end
    end


    for ku = 1:ngroups(gmap)
        for k in dataindices(gmap,ku)
            @test time[k] in rolling_classes[ku]
        end
    end
end


time = collect(DateTime(1980,1,1):Day(1):DateTime(1980,1,20))
coordinate = time

nn = 7

rolling_classes_ = [time[max(n-3,1):min(n+3,length(time))] for n = 1:length(time)]

#rolling_classes_ = [time[max(n-3,1):min(n+3,length(time))] for n = [4,8]]
#rolling_classes_ = [time[2] .+ Day.(-3:3), time[6] .+ Day.(-3:3)]


rolling_classes = [MClass(tuple(r...)) for r in rolling_classes_]

gmap = OverlappingGroupMapping(coordinate,rolling_classes)
test_mapping(gmap,coordinate,rolling_classes)



time = collect(DateTime(1980,1,1):Day(1):DateTime(1982,1,20))
class = Dates.Month.(time)
unique_class = sort(unique(class))
gmap = GroupMapping(class,unique_class)

for ku = 1:12
    @test all(Dates.Month.(time[dataindices(gmap,ku)]) .== Dates.Month(ku))
end

for k = 1:length(time)
    @test all(unique_class[collect(groupindices(gmap,k))] .== Dates.Month(time[k]))
end

# ----
TDS = MemoryDataset

time = DateTime(2000,1,1):Day(1):DateTime(2000,12,31);
data = rand(Float32.(-9:99),10,11,length(time))
varname = "data"

fname = tempname()

data3 = Array{Union{Missing,Float32},3}(undef,size(data))
data3 .= data
data3[1,1,:] .= missing
data3[1,2,1] = missing

TDS(fname,"c") do ds
    defVar(ds,"lon",1:size(data,1),("lon",))
    defVar(ds,"lat",1:size(data,2),("lat",))
    defVar(ds,"time",time,("time",))
    defVar(ds,"data",data,("lon","lat","time"))
    defVar(ds,"data3",data3,("lon","lat","time"))
end

ds = TDS(fname)

# weekly averages
nn = 7

v = ds["data"]
gv = rolling(v,:time => 7)

data_weekly_ref = similar(data);

for reduce_fun in (sum,mean,var,std,median,maximum)
    local gr
    gr = reduce_fun(gv)

    for n = 1:length(time)
        data_weekly_ref[:,:,n] = reduce_fun(data[:,:,max(n-3,1):min(n+3,length(time))],dims=3)
    end

    for indices = [
        (:,:,:),
        (1:3,2:5,:),
        (1:3,2,:),
        (2,1:3,:),
        (1:3,2:6,2:6),
        (:,:,2),
        ]

        @test gr[indices...] ≈ data_weekly_ref[indices...]
    end
end


# test data with missing values

data3 = ds["data3"][:,:,:]
data_weekly_ref = similar(data3)
for n = 1:length(time)
    data_weekly_ref[:,:,n] = mean(data3[:,:,max(n-3,1):min(n+3,length(time))],dims=3)
end

gr = mean(rolling(ds["data3"],:time => 7))[:,:,:]

@test ismissing.(gr) == ismissing.(data_weekly_ref)
@test collect(skipmissing(gr[:])) ≈ collect(skipmissing(data_weekly_ref[:]))

rolling_classes = [(t - Dates.Month(1)):Day(1):(t + Dates.Month(1)) for t in time]
extrema(length.(rolling_classes))
data_3monthly = mean(rolling(ds["data"],:time => rolling_classes))[:,:,:];


data_3monthly_ref = similar(data)
for n = 1:length(time)
    k = findall(in(rolling_classes[n]),time)
    data_3monthly_ref[:,:,n] = mean(data[:,:,k],dims=3)
end

@test data_3monthly_ref ≈ data_3monthly

using Test
using Dates
using Printf
#using NCDatasets
#using NCDatasets: NetCDFError, load!
using DataStructures
using CFTime
using CommonDataModel
using CommonDataModel:
    MemoryDataset,
    name,
    load!

sz = (4,5)
filename = tempname()
if isfile(filename)
    rm(filename)
end

#TDS = NCDataset
TDS = MemoryDataset

# The mode "c" stands for creating a new file (clobber)
TDS(filename,"c") do ds

    # define the dimension "lon" and "lat"
    ds.dim["lon"] = sz[1]
    ds.dim["lat"] = sz[2]

    v = defVar(ds,"small",Float64,("lon","lat"))
    @test parent(v) isa CommonDataModel.MemoryVariable
    @test parent(parent(v)) isa Array
#    @test_throws Union{NetCDFError,DimensionMismatch} v[:] = zeros(sz[1]+1,sz[2])
    @test_throws DimensionMismatch v[1:sz[1],1:sz[2]] = zeros(sz[1]+1,sz[2])
    @test_throws BoundsError v[sz[1]+1,1] = 1
    @test_throws BoundsError v[-1,1] = 1

    # variables
    for T in [UInt8,Int8,UInt16,Int16,UInt32,Int32,UInt64,Int64,Float32,Float64,
              Char,String]
    #for T in [String]
        local data
        data, scalar_data =
            if T == String
                [Char(i+60) * Char(j+60) for i = 1:sz[1], j = 1:sz[2]], "abcde"
            else
                [T(i+2*j) for i = 1:sz[1], j = 1:sz[2]], T(100)
            end

        v = defVar(ds,"var-$T",T,("lon","lat"))
        v[:,:] = data
        @test v[:,:] == data[:,:]

        # issue NCDatasets #33
        @test Array(v) == data

        @test v[2,:] == data[2,:]

        @test v[:,3] == data[:,3]

        @test v[2,3] == data[2,3]

        # ignore extra index
        @test v[2,3,1,1] == data[2,3,1,1]

        # ignore extra index
        @test v[2:3,3,1,1] == data[2:3,3,1,1]

        @test v[[1,2,3],:] == data[[1,2,3],:]

        # cartesian indices
        ci = CartesianIndex(1,1):CartesianIndex(3,2)
        @test v.var[ci] == data[ci]

        # write scalar
        v.var[1,1] = scalar_data
        v.var[:,:] .= scalar_data
        @test all(v.var[:,:][:] .== scalar_data)

        # stridded write and read
        v[1:2:end,1:2:end] = data[1:2:end,1:2:end]
        @test all(v[1:2:end,1:2:end] .== data[1:2:end,1:2:end])
    end
end

# quick interface
TDS(filename,"c") do ds
    data = Int32[i+3*j for i = 1:sz[1], j = 1:sz[2]]
    defVar(ds,"temp",data,("lon","lat"), attrib = [
        "units" => "degree_Celsius",
        "long_name" => "Temperature"
    ])
    @test ds["temp"][:] == data[:]
    @test eltype(ds["temp"].var) == Int32
    @test ds.dim["lon"] == sz[1]
    @test ds.dim["lat"] == sz[2]

    # load in-place of Variable
    data2 = similar(data)
    load!(ds["temp"].var,data2,:,:)
    @test data2 == data

    data2 = zeros(eltype(data),sz[1],2)
    load!(ds["temp"].var,data2,:,1:2)
    @test data2 == data[:,1:2]

    data2 = zeros(eltype(data),sz[1],1)
    load!(ds["temp"].var,data2,:,1)
    @test data2[:] == data[:,1]

    # load in-place of CFVariable
    ncv = ds["temp"]
    data2 = similar(data)
    buffer = zeros(eltype(ncv.var),size(ncv));
    load!(ncv,data2,buffer,:,:)
    @test data2 == data

    # test Union{Missing,T}
    data = [missing,1.,2.]
    defVar(ds,"foo",data,("dim",), fillvalue = -9999.)
    @test fillvalue(ds["foo"]) == -9999.
    @test isequal(ds["foo"][:], data)

    # load in-place of CFVariable with fill value
    ncv = ds["foo"]
    data2 = zeros(eltype(ncv),size(ncv))
    buffer = zeros(eltype(ncv.var),size(ncv));
    load!(ncv,data2,buffer,:)
    @test isequal(data2,data)


    # test Union{Missing,T} and default fill value (issue NCDatasets #38)
    defVar(ds,"foo_default_fill_value",[missing,1.,2.],("dim",))
    @test fillvalue(ds["foo_default_fill_value"]) == fillvalue(Float64)
    @test isequal(ds["foo_default_fill_value"][:], [missing,1.,2.])


    for DT in [DateTime,
               DateTimeStandard,
               DateTimeJulian,
               DateTimeProlepticGregorian,
               DateTimeAllLeap,
               DateTimeNoLeap,
               DateTime360Day
               ]

        # test DateTime et al., array
        data_dt = [DT(2000,1,1),DT(2000,1,2),DT(2000,1,3)]
        defVar(ds,"foo_$(DT)",data_dt,("dim",))
        data_dt2 = ds["foo_$(DT)"][:]
        @test isequal(convert.(DT,data_dt2), data_dt)

        # test DateTime et al. with missing array
        data_dt = [missing,DT(2000,1,2),DT(2000,1,3)]
        defVar(ds,"foo_$(DT)_with_fill_value",data_dt,("dim",))

        data_dt2 = ds["foo_$(DT)_with_fill_value"][:]
        @test ismissing(data_dt2[1])

        @test isequal(convert.(DT,data_dt2[2:end]), data_dt[2:end])
    end

    defVar(ds,"scalar",123.)
    @test ds["scalar"][1] == 123.

    # test indexing with symbols #101
    @test ds[:scalar][1] == 123.
end

# check bounds error
filename = tempname()
TDS(filename,"c") do ds
    defVar(ds,"temp",randn(10,11),("lon","lat"))
    @test_throws ErrorException defVar(ds,"salt",randn(10,12),("lon","lat"))
end

# check error for unknown variable
filename = tempname()
TDS(filename,"c") do ds
    @test_throws KeyError ds["does_not_exist"]
end

# issue NCDatasets 23
# return type using CartesianIndex

filename = tempname()
ds = TDS(filename, "c");
ds.dim["lon"] = 5;
ds.dim["lat"] = 10;
ds.dim["time"] = Inf;

ncvar = defVar(ds, "var", Int64, ("lon", "lat", "time"));

nt = 25;
data = reshape(1:5*10*nt, 5, 10, nt);
ncvar[:,:,1:nt] = data;
close(ds);

ds = TDS(filename);
start = 1;
all(data[CartesianIndex(1, 1), start:end] .== ds["var"][CartesianIndex(1, 1), start:end])
data11 = ds["var"][CartesianIndex(1, 1), start:end]
close(ds)

@test typeof(data11[1]) == Int64

# issue NCDatasets #36

x, y = collect(1:10), collect(10:18)

filename = tempname()
TDS(filename, "c") do ds
      defDim(ds, "x", length(x))
      defVar(ds, "x", x, ("x",))
      defDim(ds, "y", length(y))
      defVar(ds, "y", y, ("y",))
end

# issue NCDatasets 155

filename = tempname()
x = 1.:0.1:10.
ds = TDS(filename,"c");
defDim(ds, "x", length(x))
ncv = defVar(ds, "x", Float64, ("x",))
ncv[:] = x
ds.attrib["x_range"] = x
close(ds)

# issue NCDatasets 180
filename = tempname()
ds = TDS(filename, "c")

sample_data = [UInt8(1),Int64(2),Float64(3.),"string",'a']
sample_data = [UInt8(1),"string"]

for data = sample_data
    local ncv, T
    T = typeof(data)

#=    ncv = defVar(ds, "$(T)_scalar1", T, ())
    ncv[] = data
    @test ncv[] == data
=#
    ncv = defVar(ds, "$(T)_scalar2", data, ())
    @test ncv[] == data
#=
    ncv = defVar(ds, "$(T)_scalar3", data)
    @test ncv[] == data=#
end
close(ds)

# issue NCDatasets 207
filename_src = tempname()
ds_src = TDS(filename_src, "c")
data = [DateTime(2000,1,1),DateTime(2000,1,2)]
v_src = defVar(ds_src,"time",data,("time",), attrib = OrderedDict(
    "units" => "days since 2000-01-01",
))

filename_dest = tempname()
ds_dest = TDS(filename_dest, "c")
v_dest = defVar(ds_dest,v_src)
v_dest[:] = v_dest[:] .+ Dates.Minute(30)
@test name(v_dest) == name(v_src)
@test v_dest[:] == v_src[:] .+ Dates.Minute(30)

close(ds_src)
close(ds_dest)

# issue NCDatasets 209
filename_src = tempname()
ds = TDS(filename_src, "c")
data = [1,2,3]
ncv = defVar(ds,"data",data,("data",))
@test isempty(ncv[Int[]])
close(ds)

# issue NCDatasets 211
filename = tempname()
ds = TDS(filename, "c")
data = [1,2,3]
ncv = defVar(ds,"data",data,("data",))
data2 = zeros(Int,1)
# data2 too small
@test_throws DimensionMismatch load!(ds["data"].var,data2,:)

data2 = zeros(Int,10)
# asking too many elements
@test_throws BoundsError load!(ds["data"].var,data2,1:10)

# issue 22
filename = tempname()
ds = TDS(filename, "c")
v = defVar(ds, "a", Int32, ())
@test ndims(Array(v)) == ndims(v) == 0

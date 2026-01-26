#using NCDatasets
using Test
import CommonDataModel as CDM
using DataStructures
using Dates
import CommonDataModel:
    AbstractDataset,
    AbstractVariable,
    Attributes,
    Dimensions,
    defDim,
    defVar,
    nonuniontype,
    fillvalue

#include("memory_dataset.jl")

TDS = MemoryDataset

fname = tempname()
ds = MemoryDataset(fname,"c")

data = Array{Union{Missing,Float32},2}(undef,10,10)
data .= 3
data[2,2] = missing

add_offset = 12.
scale_factor = 45
fill_value = 9999.f0

v = defVar(ds,"temp",data,("lon","lat"),attrib = Dict(
    "_FillValue" => fill_value,
    "add_offset" => add_offset,
    "scale_factor" => scale_factor))

v.var[1,1] = 1

@test v[1,1] ≈ scale_factor * v.var[1,1] + add_offset
@test ismissing(v[2,2])
@test v.attrib["_FillValue"] == fill_value
@test fillvalue(v) == fill_value

@test collect(CDM.dimnames(v)) == ["lon","lat"]
@test CDM.dim(v,"lon") == 10

io = IOBuffer()
CDM.show(io,"text/plain",v)
@test occursin("Attributes",String(take!(io)))

v = @test_logs (:warn,r"numeric") defVar(ds,"temp2",data,("lon","lat"),attrib = Dict("missing_value" => "bad_idea"))

for sample_data = ( -100:100,
                    'a':'z',
                    ["AB","CD","EF"],
                    [NaN; 1:10],
                    )

    local io, data, fill_value, mv, md, add_offset, scale_factor
    local data2

    fill_value = sample_data[1]
    data = rand(sample_data[2:end],3,4)

    md = MemoryDataset()
    CDM.defDim(md,"lon",size(data,1))
    CDM.defDim(md,"lat",size(data,2))
    add_offset = 1
    scale_factor = 10

    mv = CDM.defVar(md,"data",eltype(data),("lon","lat"), attrib =
        OrderedDict(
            "_FillValue" => fill_value,
            "add_offset" => add_offset,
            "scale_factor" => scale_factor,
        ))

    mv.var[:,:] .= data

    @test "lon" in CDM.dimnames(mv)
    @test CDM.name(mv) == "data"

    md["data"][1,1] = missing
    @test ismissing(md["data"][1,1])
    @test md["data"].var[1,1] === fill_value

    md["data"][1:2,1:2] .= missing
    @test all(ismissing.(md["data"][1:2,1:2]))
    @test all(md["data"].var[1:2,1:2] .=== fill_value)

    if eltype(data) <: Number
        mv.var[3,3] = 3
        @test mv[3,3] == scale_factor * 3 + add_offset

        mv[3,3] = scale_factor * 4 + add_offset
        @test mv.var[3,3] == 4
    elseif eltype(data) == Char
        # ignore scale_factor and add_offset
        mv.var[3,3] = 'z'
        @test mv[3,3] == 'z'

        mv[3,3] = 'y'
        @test mv.var[3,3] == 'y'
    end

    # defVar(ds,name,data,dimnames)

    data2 = replace(data,fill_value => missing)
    mv = CDM.defVar(md,"data2",data2,("lon","lat"), attrib =
        OrderedDict(
            "_FillValue" => fill_value
        ))

    @test all(mv[:,:] .=== data2)
end


# time

sample_data =  -100:100
data = rand(sample_data,3,4)

md = MemoryDataset()
CDM.defDim(md,"lon",size(data,1))
CDM.defDim(md,"lat",size(data,2))


mv = CDM.defVar(md,"data",eltype(data),("lon","lat"), attrib = OrderedDict{String,Any}(
    "units" => "days since 2000-01-01"))

CDM.defAttrib(mv,"foo","bar")

@test CDM.attrib(mv,"foo") == "bar"

mv.var[:,:] .= data

@test CDM.dim(md,"lon") == size(data,1)
@test CDM.dim(mv,"lon") == size(data,1)

CDM.defAttrib(md,"history", "lala")

@test "lon" in CDM.dimnames(mv)
@test CDM.name(mv) == "data"

time_origin = DateTime(2000,1,1)
@test md["data"][1,1] == time_origin + Dates.Millisecond(data[1,1]*24*60*60*1000)


md["data"][1,2] = DateTime(2000,2,1)
@test md["data"].var[1,2] == Dates.value(md["data"][1,2] - time_origin) ÷ (24*60*60*1000)


@test CDM.dataset(md["data"]) == md

md.attrib["history"] == "lala"

@test haskey(md.attrib,"history")

@test get(md.attrib,"foooo","bar") == "bar"

@test collect(keys(md.attrib)) == ["history"]

io = IOBuffer()
show(io,md.attrib)
str = String(take!(io))
@test occursin("history",str)

show(io,md.dim)
str = String(take!(io))
@test occursin("lon",str)

@test_logs (:warn, r".*analysis.*")    CDM.defVar(md,"data2",eltype(data),("lon","lat"), attrib = OrderedDict{String,Any}(    "units" => "days since analysis"));

close(md)

# Alternative to Missing for NetCDF fillvalue

fname = tempname()
data = randn(3,4)
fv = 9999
data[2,2] = fv

ds = TDS(fname,"c")
defDim(ds,"lon",size(data,1))
defDim(ds,"lat",size(data,2))
ncv = defVar(ds,"data",Float64,("lon","lat"),fillvalue = fv)
ncv.var[:,:] = data

ncv = CDM.cfvariable(ds,"data",maskingvalue = NaN)
@test eltype(ncv) == Float64
@test ncv[1,1] == data[1,1]
@test isnan(ncv[2,2])

ncv = CDM.cfvariable(ds,"data",maskingvalue = 42)
@test eltype(ncv) == Float64
@test ncv[1,1] == data[1,1]
@test ncv[2,2] == 42

ncv = CDM.cfvariable(ds,"data",maskingvalue = nothing)
@test eltype(ncv) == Union{Nothing,Float64}
@test ncv[1,1] == data[1,1]
@test isnothing(ncv[2,2])

# custom singelton type
struct MySingleton end
const mysingleton = MySingleton()
Base.promote_rule(::Type{MySingleton},S::Type) = Union{MySingleton,S}

ncv = CDM.cfvariable(ds,"data",maskingvalue = mysingleton)
@test eltype(ncv) == Union{MySingleton,Float64}
@test ncv[1,1] == data[1,1]
@test ncv[2,2] === mysingleton
close(ds)

add_offset = -1.0
scale_factor = 0.1

attributes = Dict("add_offset" => add_offset,"scale_factor" => scale_factor)



for Tbase in (UInt8, Int8, Float32, Float64, Int64, Char, String)
    for _maskingvalue in (NaN,NaN32,missing,42,nothing,mysingleton)
        local fname, data, fv, ds, ncv, T, varattrib
        fname = tempname()
        T = promote_type(Tbase,typeof(_maskingvalue))
        data = Matrix{T}(undef,(3,4))

        if Tbase == String
            data[:,:] = rand(string.('a':'z'),size(data))
            fv = "X"
        else
            data[:,:] = rand(Tbase(0):Tbase(10),size(data))
            fv = Tbase(99)
        end

        data[2,2] = _maskingvalue
        ds = TDS(fname,"c", maskingvalue = _maskingvalue)

        varattrib =
            if Tbase <: Number
                attributes
            else
                []
            end

        ncv = defVar(ds,"data",data,("lon","lat"),fillvalue = fv, attrib = varattrib)
        @test CDM.maskingvalue(ds) === _maskingvalue
        @test ncv.var[2,2] == fv
        @test ncv[2,2] === _maskingvalue
        if Tbase <: Number
            @test ncv.var[1,1] * scale_factor + add_offset ≈ data[1,1]
            @test ncv[1,1] ≈ data[1,1]
        else
            @test ncv[1,1] == data[1,1]
        end
    end
end

@test nonuniontype(Nothing,Union{Float64,Nothing}) == Float64
@test nonuniontype(Nothing,Union{Nothing,Float64}) == Float64

@test nonuniontype(Missing,Union{Float64,Missing}) == Float64
@test nonuniontype(Missing,Union{Missing,Float64}) == Float64



fname = tempname()
data = rand(4,5)
ds = TDS(fname,"c", maskingvalue = NaN)
ncv = defVar(ds,"data",data,("lon","lat"))
@test !haskey(ncv.attrib,"_FillValue")

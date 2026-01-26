using DataStructures
using Test
sz = (4,5)
filename = tempname()
#filename = "/tmp/mytest.nc"

#TDS = NCDatasets.NCDataset
TDS = MemoryDataset

TDS(filename,"c") do ds

    ds.dim["lon"] = sz[1]
    ds.dim["lat"] = sz[2]

    v = defVar(ds,"temperature",Float32,("lon","lat"),
                          attrib = ["long_name" => "Temperature",
                                    "test_vector_attrib" => [1,2,3]])

    # write attributes
    v.attrib["units"] = "degree Celsius"
    v.attrib["comment"] = "this is a string attribute with unicode Ω ∈ ∑ ∫ f(x) dx "

    # check presence of attribute
    @test haskey(v.attrib,"comment")
    @test v.attrib["long_name"] == "Temperature"
    @test v.attrib[:long_name] == "Temperature"
    @test v.attrib["test_vector_attrib"] == [1,2,3]
    @test v.attrib["comment"] == "this is a string attribute with unicode Ω ∈ ∑ ∫ f(x) dx "

    @test get(v.attrib,"does-not-exist","default") == "default"
    @test get(v.attrib,"units","default") == "degree Celsius"

    # test deletion of attributes
    v.attrib["todelete"] = "foobar"
    @test haskey(v.attrib,"todelete")
    delete!(v.attrib,"todelete")
    @test !haskey(v.attrib,"todelete")

    for T in [UInt8,Int8,UInt16,Int16,UInt32,Int32,UInt64,Int64,Float32,Float64,
              String,Char]
        # scalar attribute
        name = "scalar-attrib-$T"

        refdata =
            if T == Char
                'a'
            elseif T == String
                "abc"
            else
                123
            end

        v.attrib[name] = T(refdata)
        attval = v.attrib[name]

        @test typeof(attval) == T
        @test attval == refdata

        # vector attribute
        name = "vector-attrib-$T"

        refvecdata =
            if T == Char
                ['a','b']
            elseif T == String
                ["abc","xyz"]
            else
                [1,2,3,4]
            end


        attval = T.(refvecdata)
        attval = attval

        @test eltype(attval) == T
        @test attval == refvecdata
    end

    # symbols in the attrib dict
    foo = defVar(ds,"foovar",Int64,("lon","lat"),
                            attrib = [:long_name => "foo variable"])
    @test foo.attrib["long_name"] == "foo variable"

end

filename = tempname()

TDS(filename,"c") do ds
    # test deletion of attributes
    ds.attrib["todelete"] = "foobar"
end

TDS(filename,"a") do ds
    @test haskey(ds.attrib,"todelete")
    delete!(ds.attrib,"todelete")
    @test !haskey(ds.attrib,"todelete")
end


# NCDatasets issue #241

filename = tempname()
ds = TDS(filename,"c")
data = randn(4,5)

defVar(ds,"temperature",data,("lon","lat"),
       attrib = OrderedDict(:long_name => "Temperature",
                 :test_vector_attrib => [1,2,3]))

defVar(ds,:temperature2,data,(:lon,:lat),
       attrib = OrderedDict(:long_name => "Temperature",
                 :test_vector_attrib => [1,2,3]))

defVar(ds,:temperature3,data,(:lon,:lat),
       attrib = Dict(:long_name => "Temperature",
                 :test_vector_attrib => [1,2,3]))

defVar(ds,:temperature4,data,(:lon,:lat),
       attrib = (long_name = "Temperature",
                 test_vector_attrib = [1,2,3]))

for (vn,var) in ds
    @test startswith(vn,"temperature")
    @test var.attrib["long_name"] == "Temperature"
    @test var.attrib["test_vector_attrib"] == [1,2,3]
end

"""
    CommonDatamodel.attribnames(ds::Union{AbstractDataset,AbstractVariable})

Return an iterable of all attribute names in `ds`.
"""
attribnames(ds::Union{AbstractDataset,AbstractVariable}) = ()


"""
    CommonDatamodel.attrib(ds::Union{AbstractDataset,AbstractVariable},attribname::SymbolOrString)

Return the value of the attribute `attribname` in the data set `ds`.
"""
function attrib(ds::Union{AbstractDataset,AbstractVariable},attribname::SymbolOrString)
    error("no attributes $attribname in $(path(ds))")
end

"""
    CommonDatamodel.defAttrib(ds::Union{AbstractDataset,AbstractVariable},name::SymbolOrString,data)

Create an attribute with the name `attrib` in the data set or variable `ds`.
"""
function defAttrib(ds::AbstractDataset,name::SymbolOrString,data)
    error("unimplemnted for abstract type")
end


"""
    CommonDatamodel.delAttrib(ds::Union{AbstractDataset,AbstractVariable},name::SymbolOrString,data)

Deletes an attribute with the name `attrib` in the data set or variable `ds`.
"""
function delAttrib(ds::Union{AbstractDataset,AbstractVariable},name::SymbolOrString,data)
    error("unimplemnted for abstract type")
end


attribs(ds::Union{AbstractDataset,AbstractVariable}) =
    OrderedDict((dn,attrib(ds,dn)) for dn in attribnames(ds))



#=
    CommonDatamodel.show_attrib(io,a)

Print a list all attributes (key/values pairs) in `a` to IO stream `io`.
The IO property `:level` is used for indentation.
=#
function show_attrib(io,a)
    level = get(io, :level, 0)
    indent = " " ^ level

    # need to know ds from a
    #if !isopen(ds)
    #    print(io,"Dataset attributes (file closed)")
    #    return
    #end

    try
        # use the same order of attributes than in the dataset
        for (attname,attval) in a
            print(io,indent,@sprintf("%-20s = ",attname))
            printstyled(io, @sprintf("%s",attval),color=attribute_color[])
            print(io,"\n")
        end
    catch err
        @debug "error while printing" err
        print(io,"Dataset attributes (file closed)")
    end
end


"""
    Base.keys(a::Attributes)

Return a list of the names of all attributes.
"""
Base.keys(a::Attributes) = attribnames(a.ds)


"""
    getindex(a::Attributes,name::SymbolOrString)

Return the value of the attribute called `name` from the
attribute list `a`. Generally the attributes are loaded by
indexing, for example:

```julia
using NCDatasets
ds = NCDataset("file.nc")
title = ds.attrib["title"]
```
"""
Base.getindex(a::Attributes,name) = attrib(a.ds,name)


"""
    Base.setindex!(a::Attributes,data,name::SymbolOrString)

Set the attribute called `name` to the value `data` in the
attribute list `a`. `data` can be a vector or a scalar. A scalar
is handeld as a vector with one element in the NetCDF data model.

Generally the attributes are defined by indexing, for example:

```julia
ds = NCDataset("file.nc","c")
ds.attrib["title"] = "my title"
close(ds)
```
"""
Base.setindex!(a::Attributes,data,name) = defAttrib(a.ds,name,data)

Base.show(io::IO,a::Attributes) = show_attrib(io,a)



"""
    Base.delete!(a::Attributes, name)

Delete the attribute `name` from the attribute list `a`.
"""
Base.delete!(a::Attributes,name::SymbolOrString) = delAttrib(a.ds,name)

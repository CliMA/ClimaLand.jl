
"""
    CommonDatamodel.path(ds::AbstractDataset)

File path of the data set `ds`.
"""
path(ds::AbstractDataset) = ""


Base.close(ds::AbstractDataset) = nothing

"""
    CommonDatamodel.name(ds::AbstractDataset)

Name of the group of the data set `ds`. For a data set containing
only a single group, this will be always the root group `"/"`.
"""
name(ds::AbstractDataset) = "/"

"""
    pds = CommonDatamodel.parentdataset(ds::AbstractDataset)

The data set `pds` containing `ds` as a sub-group.
`pds` is nothing for the root group.
"""
parentdataset(ds::AbstractDataset) = nothing

sync(ds::AbstractDataset) = nothing

"""
    CommonDatamodel.groupnames(ds::AbstractDataset)

All the subgroup names of the data set `ds`. For a data set containing
only a single group, this will be an empty vector of `String`.
"""
groupnames(ds::AbstractDataset) = []


"""
    CommonDatamodel.group(ds::AbstractDataset,groupname::SymbolOrString)

Return the sub-group data set with the name `groupname`.
"""
function group(ds::AbstractDataset,groupname::SymbolOrString)
    error("no group $groupname in $(path(ds))")
end

"""
    group = CommonDatamodel.defGroup(ds::AbstractDataset,name::SymbolOrString)

Create an empty sub-group with the name `name` in the data set `ds`.
The `group` is a sub-type of `AbstractDataset`.
"""
function defGroup(ds::AbstractDataset,name::SymbolOrString)
    error("unimplemented for abstract type")
end

"""
    CommonDatamodel.groups(ds::AbstractDataset)

Return all sub-group data as a dict-like object.
"""
groups(ds::AbstractDataset) =
    OrderedDict((dn,group(ds,dn)) for dn in groupnames(ds))


"""
    CommonDatamodel.unlimited(ds::AbstractDataset)

Iterator of strings with the name of the unlimited dimension.
"""
unlimited(ad::AbstractDataset) = ()

Base.isopen(ds::AbstractDataset) = true

iswritable(ds::AbstractDataset) = false

function Base.show(io::IO,ds::AbstractDataset)
    level = get(io, :level, 0)
    indent = " " ^ level

    if !isopen(ds)
        print(io,"closed Dataset")
        return
    end

    dspath = path(ds)
    printstyled(io, indent, "Dataset: ",dspath,"\n", color=section_color[])

    print(io,indent,"Group: ",name(ds),"\n")
    print(io,"\n")

    # show dimensions
    if length(dimnames(ds)) > 0
        show_dim(io, dims(ds))
        print(io,"\n")
    end

    varnames = keys(ds)

    if length(varnames) > 0
        printstyled(io, indent, "Variables\n",color=section_color[])

        for name in varnames
            show(IOContext(io,:level=>level+2),ds[name])
            print(io,"\n")
        end
    end

    # global attribues
    if length(attribnames(ds)) > 0
        printstyled(io, indent, "Global attributes\n",color=section_color[])
        show_attrib(IOContext(io,:level=>level+2),attribs(ds));
    end

    # groups
    gnames = groupnames(ds)

    if length(gnames) > 0
        printstyled(io, indent, "Groups\n",color=section_color[])
        for groupname in gnames
            show(IOContext(io,:level=>level+2),group(ds,groupname))
        end
    end

end



maskingvalue(ds::AbstractDataset) = missing

"""
    v = getindex(ds::AbstractDataset, varname::SymbolOrString)

Return the variable `varname` in the dataset `ds` as a
`CFVariable`. The following CF convention are honored when the
variable is indexed:
* `_FillValue` or `missing_value` (which can be a list) will be returned as `missing`.
* `scale_factor` and `add_offset` are applied (output = `scale_factor` * `data_in_file` +  `add_offset`)
* time variables (recognized by the units attribute and possibly the calendar attribute) are returned usually as
  `DateTime` object. Note that `CFTime.DateTimeAllLeap`, `CFTime.DateTimeNoLeap` and
  `CF.TimeDateTime360Day` cannot be converted to the proleptic gregorian calendar used in
  julia and are returned as such. (See [`CFTime.jl`](https://github.com/JuliaGeo/CFTime.jl)
  for more information about those date types.) If a calendar is defined but not among the
  ones specified in the CF convention, then the data in the file is not
  converted into a date structure.

A call `getindex(ds, varname)` is usually written as `ds[varname]`.

If variable represents a cell boundary, the attributes `calendar` and `units` of the related variables are used, if they are not specified. For example:

```
dimensions:
  time = UNLIMITED; // (5 currently)
  nv = 2;
variables:
  double time(time);
    time:long_name = "time";
    time:units = "hours since 1998-04-019 06:00:00";
    time:bounds = "time_bnds";
  double time_bnds(time,nv);
```

In this case, the variable `time_bnds` uses the units and calendar of `time`
because both variables are related thought the bounds attribute following the CF conventions.

See also [`cfvariable(ds, varname)`](@ref).
"""
function Base.getindex(ds::AbstractDataset,varname::SymbolOrString)
    return cfvariable(ds, varname)
end


function Base.setindex!(ds::AbstractDataset,data::AbstractVariable,varname::SymbolOrString)
    return defVar(ds, varname, data)
end

function Base.keys(ds::AbstractDataset)
    return varnames(ds)
end

function Base.haskey(ds::AbstractDataset,varname)
    return Symbol(varname) in Symbol.(keys(ds))
end

"""
    varbyattrib(ds, attname = attval)

Returns a list of variable(s) which has the attribute `attname` matching the value `attval`
in the dataset `ds`.
The list is empty if the none of the variables has the match.
The output is a list of `CFVariable`s.

# Examples

Load all the data of the first variable with standard name "longitude" from the
NetCDF file `results.nc`.

```julia-repl
julia> ds = NCDataset("results.nc", "r");
julia> data = varbyattrib(ds, standard_name = "longitude")[1][:]
```

"""
function varbyattrib(ds::Union{AbstractDataset,AbstractVariable}; kwargs...)
    # Start with an empty list of variables
    varlist = []

    # Loop on the variables
    for v in keys(ds)
        var = ds[v]

        matchall = true

        for (attsym,attval) in kwargs
            attname = String(attsym)

            # Check if the variable has the desired attribute
            if attname in attribnames(var)
                # Check if the attribute value is the selected one
                if attrib(var,attname) != attval
                    matchall = false
                    break
                end
            else
                matchall = false
                break
            end
        end

        if matchall
            push!(varlist, var)
        end
    end

    return varlist
end

"""
    var = getindex(ds::Union{AbstractDataset,AbstractVariable},cfname::CFStdName)

Return the NetCDF variable `var` with the standard name `cfname` from a
dataset. If the first argument is a variable, then the search is limited to
all variables with the same dimension names.
"""
function Base.getindex(ds::Union{AbstractDataset,AbstractVariable},n::CFStdName)
    ncvars = varbyattrib(ds, standard_name = String(n.name))
    if length(ncvars) == 1
        return ncvars[1]
    else
        throw(KeyError("variable with standard_name attribute equal to $(n.name) ($(length(ncvars)) matches)"))
    end
end

"""
    names = keys(g::Groups)

Return the names of all subgroubs of the group `g`.
"""
Base.keys(groups::Groups) = groupnames(groups.ds)

"""
    group = getindex(g::Groups,groupname::AbstractString)

Return the NetCDF `group` with the name `groupname` from the parent
group `g`.

For example:

```julia
ds = NCDataset("results.nc", "r");
forecast_group = ds.group["forecast"]
forecast_temp = forecast_group["temperature"]
```
"""
Base.getindex(groups::Groups,name) = group(groups.ds,name)


# Initialize the ds._boundsmap variable
function initboundsmap!(ds)
    empty!(ds._boundsmap)
    for vname in keys(ds)
        v = variable(ds,vname)
        bounds = get(v.attrib,"bounds",nothing)

        # see https://github.com/Alexander-Barth/NCDatasets.jl/issues/265
        if bounds isa String
            ds._boundsmap[bounds] = vname
        end
    end
end


"""
    write(dest::AbstractDataset, src::AbstractDataset; include = keys(src), exclude = [])

Write the variables of `src` dataset into an empty `dest` dataset (which must be opened in mode `"a"` or `"c"`).
The keywords `include` and `exclude` configure which variable of `src` should be included
(by default all), or which should be `excluded` (by default none).

If the first argument is a file name, then the dataset is open in create mode (`"c"`).

This function is useful when you want to save the dataset from a multi-file dataset.

To save a subset, one can use the view function `view` to virtually slice
a dataset:

## Example

```
NCDataset(fname_src) do ds
    write(fname_slice,view(ds, lon = 2:3))
end
```

All variables in the source file `fname_src` with a dimension `lon` will be sliced
along the indices `2:3` for the `lon` dimension. All attributes (and variables
without a dimension `lon`) will be copied over unmodified.
"""
function Base.write(dest::AbstractDataset, src::AbstractDataset;
                    include = keys(src),
                    exclude = String[],
                    _ignore_checksum = false,
                    )


    unlimited_dims = unlimited(src)

    for (dimname,dimlength) in dims(src)
        isunlimited = dimname in unlimited_dims

        # if haskey(dest.dim,dimname)
        #     # check length
        #     if (dest.dim[dimname] !== src.dim[dimname]) && !isunlimited
        #         throw(DimensionMismatch("length of the dimensions $dimname are inconstitent in files $(path(dest)) and $(path(src))"))
        #     end
        # else
            if isunlimited
                defDim(dest, dimname, Inf)
            else
                defDim(dest, dimname, dimlength)
            end
        # end
    end

    # loop over variables
    for varname in include
        (varname ∈ exclude) && continue
        @debug "Writing variable $varname..."

        kwargs =
            if _ignore_checksum
                (checksum = nothing,)
            else
                ()
            end

        defVar(dest,src[varname]; kwargs...)
    end

    # loop over all global attributes
    for (attribname,attribval) in attribs(src)
        dest.attrib[attribname] = attribval
    end

    # loop over all groups
    for (groupname,groupsrc) in groups(src)
        groupdest = defGroup(dest,groupname)
        write(groupdest,groupsrc)
    end

    return dest
end


@inline function Base.getproperty(ds::Union{AbstractDataset,AbstractVariable},name::Symbol)
    if (name == :attrib) && !hasfield(typeof(ds),name)
        return Attributes(ds)
    elseif (name == :dim) && !hasfield(typeof(ds),name)
        return Dimensions(ds)
    elseif (name == :group) && !hasfield(typeof(ds),name) && (ds isa AbstractDataset)
        return Groups(ds)
    else
        return getfield(ds,name)
    end
end

@inline function Base.propertynames(ds::Union{AbstractDataset,AbstractVariable},private::Bool=false)
    names = fieldnames(typeof(ds))

    if ds isa AbstractDataset
        return (names...,:attrib,:dim,:group)
    else
        return (names...,:attrib,:dim)
    end
end

for (item_color,default) in (
    (:section_color, :red),
    (:attribute_color, :cyan),
    (:variable_color, :green),
)

    item_color_str = String(item_color)
    item_str = split(item_color_str,"_")[1]
    default_str = String(default)

    @eval begin
        $item_color = Ref(Symbol(load_preference(CommonDataModel,$(item_color_str), $(QuoteNode(default)))))

        """
        CommonDataModel.set_$($item_color_str)(color::Symbol)

Set the $($item_str) color. The default color is `$($default_str)`.
"""
        function $(Symbol(:set_,item_color))(color::Symbol)
            @set_preferences!($(item_color_str) => String(color))
            $item_color[] = color
        end

    end
end

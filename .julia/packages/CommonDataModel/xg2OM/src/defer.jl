iswritable(dds::DeferDataset) = dds.r.mode != "r"

function metadata(ds::AbstractDataset)
    # dimensions

    dim = OrderedDict()
    unlimited_dims = unlimited(ds)

    for (dimname,dimlen) in ds.dim
        dim[dimname] = Dict(
            :name => dimname,
            :length => dimlen,
            :unlimited => dimname in unlimited_dims
        )
    end

    # variables
    vars = OrderedDict()

    for (varname,ncvar) in ds
        storage,chunksizes = chunking(ncvar.var)
        isshuffled,isdeflated,deflatelevel = deflate(ncvar.var)
        checksummethod = checksum(ncvar.var)

        vars[varname] = OrderedDict(
            :name => varname,
            :size => size(ncvar),
            :eltype => eltype(ncvar.var),
            :attrib => OrderedDict(ncvar.attrib),
            :dimensions => dimnames(ncvar),
            :chunksize => chunksizes,
            :storage => storage,
            #:fillvalue => fillvalue(ncvar.var),
            :shuffle => isshuffled,
            :deflate => isdeflated,
            :deflatelevel => deflatelevel,
            :checksummethod => checksummethod,
        )
    end

    group = OrderedDict()
    for (groupname,ncgroup) in ds.group
        group[groupname] = metadata(ncgroup)
    end

    return OrderedDict(
        :dim => dim,
        :var => vars,
        :attrib => OrderedDict(ds.attrib),
        :group => group
        )
end


function DeferDataset(TDS,r::Resource,groupname::String,data::OrderedDict)
   _boundsmap = Dict{String,String}()
   dds = DeferDataset{TDS}(r,groupname,data,_boundsmap)
   if (r.mode == "r")
       initboundsmap!(dds)
   end
   return dds
end

function DeferDataset(TDS,filename::AbstractString,mode::AbstractString,info::OrderedDict; kwargs...)
    r = Resource(filename,mode,Dict(kwargs),info)
    groupname = "/"
    return DeferDataset(TDS,r,groupname,info)
end

function DeferDataset(TDS,filename::AbstractString,mode = "r"; kwargs...)
    TDS(filename,mode) do ds
        info = metadata(ds)
        r = Resource(filename,mode,Dict(kwargs),info)
        groupname = "/"
        return DeferDataset(TDS,r,groupname,info)
    end
end


# files are not suppose to be open where using DeferDataset
close(dds::DeferDataset) = nothing
groupname(dds::DeferDataset) = dds.groupname
path(dds::DeferDataset) = dds.r.filename
varnames(dds::DeferDataset) = collect(keys(dds.data[:var]))

function maskingvalue(dds::DeferDataset{TDS}) where TDS
    TDS(dds.r.filename,dds.r.mode; dds.r.args...) do ds
        maskingvalue(ds)
    end
end

function Variable(f::Function, dv::DeferVariable{T,N,TDS}) where {T,N,TDS}
    TDS(dv.r.filename,dv.r.mode; dv.r.args...) do ds
        f(variable(ds,dv.varname))
    end
end

function variable(dds::DeferDataset{TDS},varname::AbstractString) where TDS
    data = get(dds.data[:var],varname,nothing)
    if data == nothing
        error("Dataset $(dds.r.filename) does not contain the variable $varname")
    end
    T = data[:eltype]
    N = length(data[:dimensions])

    return DeferVariable{T,N,TDS}(dds.r,varname,data)
end

variable(dds::DeferDataset,varname::Symbol) = variable(dds,string(varname))

dataset(dv::DeferVariable{T,N,TDS}) where {T,N,TDS} =
    DeferDataset(TDS,dv.r.filename,dv.r.mode; dv.r.args...)

function Base.getindex(dv::DeferVariable,indexes::Union{Int,Colon,AbstractRange{<:Integer}}...)
    Variable(dv) do v
        return v[indexes...]
    end
end


Base.size(dv::DeferVariable) = dv.data[:size]
dimnames(dv::DeferVariable) = dv.data[:dimensions]
name(dv::DeferVariable) = dv.varname

#----------------------------------------------

dimnames(dds::DeferDataset) = collect(keys(dds.r.metadata[:dim]))

dim(dds::DeferDataset,name::SymbolOrString) =
    dds.r.metadata[:dim][String(name)][:length]

unlimited(dd::DeferDataset) = [dimname for (dimname,dim) in dd.data[:dim] if dim[:unlimited]]


attribnames(dds::DeferDataset) = collect(keys(dds.r.metadata[:attrib]))
attrib(dds::DeferDataset,name::SymbolOrString) = dds.r.metadata[:attrib][String(name)]


attribnames(dv::DeferVariable) = collect(keys(dv.data[:attrib]))
attrib(dv::DeferVariable,name::SymbolOrString) = dv.data[:attrib][String(name)]

#------------------------------------------------

groupnames(dds::DeferDataset) = collect(keys(dds.data[:group]))

function group(dds::DeferDataset{TDS},name::SymbolOrString) where TDS
    data = dds.data[:group][String(name)]
    return DeferDataset(TDS,dds.r,String(name),data)
end


_storage_attributes(dv) = dv.r.metadata[:var][name(dv)]

function chunking(dv::DeferVariable)
    sa = _storage_attributes(dv)
    return sa[:storage],sa[:chunksize]
end

function deflate(dv::DeferVariable)
    sa = _storage_attributes(dv)
    return sa[:shuffle],sa[:deflate],sa[:deflatelevel]
end

checksum(dv::DeferVariable) = _storage_attributes(dv)[:checksummethod]

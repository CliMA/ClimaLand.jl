import CommonDataModel as CDM

Base.size(RA::ResizableArray) = size(RA.A)
Base.getindex(RA::ResizableArray,inds...) = getindex(RA.A,inds...)

@inline _checkbounds(RA::ResizableArray,inds...) = all(minimum.(inds) .> 0)

@inline Base.checkbounds(::Type{Bool},RA::ResizableArray,inds::Union{Integer,AbstractVector{<:Integer}}...) =
    _checkbounds(RA,inds...)

@inline checkbounds(::Type{Bool}, itp::ResizableArray, inds::LogicalIndex) =
    _checkbounds(RA,inds...)

@inline checkbounds(::Type{Bool}, itp::ResizableArray, inds::LogicalIndex{<:Any,<:AbstractVector{Bool}})  =
    _checkbounds(RA,inds...)

# ambiguity julia 1.11
@inline checkbounds(::Type{Bool}, itp::ResizableArray, inds::AbstractVector{Bool}) =
    _checkbounds(RA,inds...)


function grow!(RA::ResizableArray{T,N},new_size) where {T,N}
    # grow
    oldA = RA.A
    RA.A = Array{T,N}(undef,new_size)
    RA.A .= RA.fillvalue
    RA.A[axes(oldA)...] = oldA
end

function Base.setindex!(RA::ResizableArray{T,N}, value, inds::Vararg{Int, N}) where {T,N}
    sz = max.(size(RA),inds)
    if sz != size(RA)
        grow!(RA,sz)
    end
    RA.A[inds...] = value
end


function _root(ds::Union{MemoryVariable,MemoryDataset})
    if isnothing(ds.parent_dataset)
        return ds
    else
        return _root(ds.parent_dataset)
    end
end

function grow_unlimited_dimension(ds,dname,len)
    if haskey(ds.dimensions,dname)
        ds.dimensions[dname] = len
    end

    for (varname,var) in ds.variables
        new_size = ntuple(ndims(var)) do j
            if dimnames(var)[j] == dname
                len
            else
                size(var,j)
            end
        end

        if new_size != size(var)
            grow!(var.data,new_size)
        end
    end

    for (groupname,group) in ds._group
        grow_unlimited_dimension(group,dname,len)
    end
end

Base.parent(v::MemoryVariable) = v.data
Base.size(v::MemoryVariable) = size(parent(v))
Base.getindex(v::MemoryVariable,ij::TIndices...) = parent(v)[ij...]
function Base.setindex!(v::MemoryVariable,data,ij...)
    sz = size(v)
    parent(v)[ij...] = data

    root = _root(v)
    for idim = findall(size(v) .> sz)
        dname = v.dimnames[idim]
        grow_unlimited_dimension(v.parent_dataset,dname,size(v,idim))
    end
    return data
end

CDM.load!(v::MemoryVariable,buffer,ij...) = buffer .= view(parent(v),ij...)
CDM.name(v::Union{MemoryVariable,MemoryDataset}) = v.name
CDM.dimnames(v::MemoryVariable) = v.dimnames
CDM.dataset(v::MemoryVariable) = v.parent_dataset

CDM.variable(md::MemoryDataset,varname::SymbolOrString) = md.variables[String(varname)]
CDM.dimnames(md::MemoryDataset) = keys(md.dimensions)
CDM.maskingvalue(md::MemoryDataset) = md.maskingvalue

function CDM.unlimited(md::MemoryDataset)
    ul = md.unlimited
    if md.parent_dataset != nothing
        append!(ul,unlimited(md.parent_dataset))
    end
    return ul
end

function _dim(md::MemoryDataset,name::SymbolOrString)
    if haskey(md.dimensions,String(name))
        return md.dimensions[String(name)]
    elseif md.parent_dataset !== nothing
        return _dim(md.parent_dataset,name)
    end
    return nothing
end

function CDM.dim(md::MemoryDataset,name::SymbolOrString)
    len = _dim(md,name)
    if !isnothing(len)
        return len
    else
        error("dimension $name not found")
    end
end

CDM.varnames(ds::MemoryDataset) = collect(keys(ds.variables))
CDM.attribnames(md::Union{MemoryDataset,MemoryVariable}) = keys(md._attrib)
CDM.attrib(md::Union{MemoryDataset,MemoryVariable},name::SymbolOrString) = md._attrib[String(name)]


CDM.groupnames(md::MemoryDataset) = collect(keys(md._group))
CDM.group(md::MemoryDataset,name::SymbolOrString) = md._group[String(name)]

function CDM.defDim(md::MemoryDataset,name::SymbolOrString,len)
    if isinf(len)
        md.dimensions[String(name)] = 0
        push!(md.unlimited,String(name))
    else
        md.dimensions[String(name)] = len
    end
end

function CDM.defVar(md::MemoryDataset,name::SymbolOrString,T::DataType,dimnames;
                    fillvalue = nothing,
                    attrib = OrderedDict{SymbolOrString,Any}(),
                    )

    sz = ntuple(i -> CDM.dim(md,dimnames[i]),length(dimnames))

    if length(intersect(dimnames,CDM.unlimited(md))) == 0
        data = Array{T,length(dimnames)}(undef,sz...)
    else
        fv =
            if !isnothing(fillvalue)
                T(fillvalue)
            elseif haskey(attrib,"_FillValue")
                T(attrib["_FillValue"])
            else
                T(0)
            end

        data_ = Array{T,length(dimnames)}(undef,sz...)
        data = ResizableArray(data_,fv)
    end

    attrib_ = OrderedDict{String,Any}()
    for (k,v) in attrib
        attrib_[String(k)] = v
    end
    mv = MemoryVariable(md,String(name),(String.(dimnames)...,), data, attrib_)

    if fillvalue !== nothing
        mv.attrib["_FillValue"] = fillvalue
    end
    md.variables[String(name)] = mv

    cfvar = md[String(name)]

    return cfvar
end

function CDM.defAttrib(md::Union{MemoryVariable,MemoryDataset},name::SymbolOrString,data)
    md._attrib[String(name)] = data
end

function CDM.delAttrib(md::Union{MemoryVariable,MemoryDataset},name::SymbolOrString)
    delete!(md._attrib,String(name))
end

function CDM.defGroup(md::MemoryDataset,name::SymbolOrString);
    md._group[String(name)] = MemoryDataset(; parent_dataset = md, name = name)
end


CDM.parentdataset(md::MemoryDataset) = md.parent_dataset
CDM.iswritable(md::MemoryDataset) = true

function MemoryDataset(; parent_dataset = nothing, name = "/",
                       attrib = OrderedDict{String,Any}(),
                       maskingvalue = missing
                       )
    return MemoryDataset(
        parent_dataset,
        name,
        OrderedDict{String,Int}(),
        OrderedDict{String,MemoryVariable}(),
        OrderedDict{String,Any}(attrib),
        String[],
        OrderedDict{String,Any}(),
        maskingvalue
    )
end

const MEMORY_DATASET_STORE = Dict{String,Any}()

function MemoryDataset(key,mode = "r"; kwargs...)
    if mode == "c"
        md = MemoryDataset(; kwargs...);
        MEMORY_DATASET_STORE[key] = md
        return md
    elseif mode == "r" || mode == "a"
        return MEMORY_DATASET_STORE[key]
    end
end


MemoryDataset(keys::AbstractArray{<:AbstractString,N}, args...; kwargs...) where N =
    MFDataset(MemoryDataset,keys, args...; kwargs...)


function MemoryDataset(f::Function,args...; kwargs...)
    ds = MemoryDataset(args...; kwargs...)
    try
        f(ds)
    finally
        close(ds)
    end
end

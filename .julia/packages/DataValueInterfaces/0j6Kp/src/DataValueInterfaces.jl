module DataValueInterfaces

"""
    nondatavaluetype(T)

For a type `T`, return the corresponding non-`DataValue` type, translating between
`Union{T, Missing}` and `DataValue{T}`.

For example, `nondatavaluetype(Int64)` returns `Int64`, while
`nondatavaluetype(DataValue{Int64})` returns `Union{Int64, Missing}`.

This generic function is owned by DataValueInterfaces.jl itself, which is the sole provider of the
default definition.
"""
function nondatavaluetype end

nondatavaluetype(::Type{T}) where {T} = T

"""
    datavaluetype(T)

For a type `T`, return the corresponding `DataValue` type, translating between 
`Union{T, Missing}` and `DataValue{T}`.

For example, `datavaluetype(Int64)` returns `Int64`, while
`datavaluetype(Union{Int64, Missing})` returns `DataValue{Int64}`.

This generic function is owned by DataValueInterfaces.jl itself, which is the sole provider of the
default definition.
"""
function datavaluetype end

datavaluetype(::Type{T}) where {T} = T

end # module

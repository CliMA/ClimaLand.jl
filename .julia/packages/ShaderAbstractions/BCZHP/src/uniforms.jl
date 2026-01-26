
# All the types native to ogl, wgl and vulkan shaders

const number_types = (Float32, Cint, Cuint, Cdouble)
const small_vecs = tuple((StaticVector{N, T} for T in number_types, N in (2, 3, 4))...)
const small_mats = tuple((Mat{R, C, T} for T in number_types, R in (2, 3, 4), C in (2, 3, 4))...)
const small_arrays = (small_vecs..., small_mats...)
const native_types = (number_types..., small_arrays...)

const NativeNumbers = Union{number_types...}
const SmallVecs = Union{small_vecs...}
const SmallMats = Union{small_mats...}
const SmallArrays = Union{small_arrays...}
const NativeTypes = Union{native_types...}

"""
    native_type(context, Int128)
Returns a native type for a non native type.
E.g. native_type(context, Int128) -> Cint
"""
native_type(context::AbstractContext, ::Type{T}) where T <: NativeTypes = T

native_type(context::AbstractContext, x::Type{T}) where {T <: Integer} = Cint
native_type(context::AbstractContext, x::Type{Union{Int16, Int8}}) = x

native_type(context::AbstractContext, x::Type{T}) where {T <: Unsigned} = Cuint
native_type(context::AbstractContext, x::Type{Union{UInt16, UInt8}})  = x

native_type(context::AbstractContext, x::Type{T}) where {T <: AbstractFloat} = Float32
native_type(context::AbstractContext, x::Type{Float16})               = x

native_type(context::AbstractContext, x::Type{T}) where T <: Normed = N0f32
native_type(context::AbstractContext, x::Type{N0f16}) = x
native_type(context::AbstractContext, x::Type{N0f8}) = x

native_type(context::AbstractContext, x::Type{<: StaticVector{N, T}}) where {T, N} = similar_type(x, native_type(context, T))
native_type(context::AbstractContext, x::Type{<: Mat{R, C, T}}) where {R, C, T} = similar_type(x, native_type(context, T))

map_t(f, tuple) = map_t(f, (), tuple)
map_t(f, result, ::Type{Tuple{}}) = Tuple{result...}
function map_t(f, result, T::Type{<: Tuple})
    map_t(
        f,
        (result..., f(Base.tuple_type_head(T))),
        Base.tuple_type_tail(T)
    )
end

native_type(context::AbstractContext) = x-> native_type(context, x)

function native_type(context::AbstractContext, ::Type{NamedTuple{Names, Types}}) where {Names, Types}
    return NamedTuple{Names, map_t(native_type(context), Types)}
end

"""
Curry form to easier pass this via e.g. map
"""
convert_uniform(context::AbstractContext) = x-> convert_uniform(context, x)

"""
All native types don't need conversion
"""
convert_uniform(context::AbstractContext, x::NativeTypes) = x

"""
Vector of native types, e.g `vec3 [4]`
"""
convert_uniform(context::AbstractContext, x::StaticVector{N, T}) where {N, T <: NativeTypes} = x

"""
Static Array with non native uniform type
"""
function convert_uniform(context::AbstractContext, x::StaticVector{N, T}) where {N, T}
    convert(similar_type(x, native_type(context, T)), x)
end

"""
Colors get special treatment to get them as vecs in shaders
"""
convert_uniform(context::AbstractContext, x::Colorant{T}) where T <: NativeNumbers = x

convert_uniform(context::AbstractContext, x::Colorant{T}) where T = mapc(native_type(context, T), x)

convert_uniform(context::AbstractContext, x::AbstractSampler{T, N}) where {T, N} = x


function convert_uniform(context::AbstractContext, x::AbstractVector{T}) where T
    return convert(Vector{native_type(context, T)}, x)
end

function convert_uniform(context::AbstractContext, x::NamedTuple{Names, Types}) where {Names, Types}
    return map(convert_uniform(context), x)
end

function convert_uniform(context::AbstractContext, x::Observable)
    return map(convert_uniform(context), x)
end
function convert_uniform(context::AbstractContext, x::T) where T
    all(t-> isbits(t), fieldtypes(T)) || error("All field types need to be isbits. Found: $(T) with $(fieldtypes(T))")
    return x
end


type_prefix(context::AbstractContext, x::Type{T}) where {T <: Union{FixedPoint, Float32, Float16}} = ""
type_prefix(context::AbstractContext, x::Type{T}) where {T <: Float64} = "d"
type_prefix(context::AbstractContext, x::Type{Cint}) = "i"
type_prefix(context::AbstractContext, x::Type{T}) where {T <: Union{Cuint, UInt8, UInt16}} = "u"

type_postfix(context::AbstractContext, x::Type{Float64}) = "dv"
type_postfix(context::AbstractContext, x::Type{Float32}) = "fv"
type_postfix(context::AbstractContext, x::Type{Cint})    = "iv"
type_postfix(context::AbstractContext, x::Type{Cuint})   = "uiv"


type_string(context::AbstractContext, x::T) where {T} = type_string(context, T)
type_string(context::AbstractContext, t::Type{Float32}) = "float"
type_string(context::AbstractContext, t::Type{Float64}) = "double"
type_string(context::AbstractContext, t::Type{Cuint}) = "uint"
type_string(context::AbstractContext, t::Type{Cint}) = "int"

function type_string(context::AbstractContext, t::Type{T}) where {T <: Union{StaticVector, Colorant}}
    return string(type_prefix(context, eltype(T)), "vec", length(T))
end

function type_string(context::AbstractContext, t::Type{<: AbstractSamplerBuffer{T}}) where T
    string(type_prefix(context, eltype(T)), "samplerBuffer")
end
is_arraysampler(t) = false
function type_string(context::AbstractContext, t::Type{<: AbstractSampler{T, D}}) where {T, D}
    str = string(type_prefix(context, eltype(T)), "sampler", D == 1 ? 2 : D, "D")
    is_arraysampler(t) && (str *= "Array")
    return str
end

function type_string(context::AbstractContext, t::Type{<: Mat})
    M, N = size(t)
    string(type_prefix(context, eltype(t)), "mat", M == N ? M : string(M, "x", N))
end
type_string(context::AbstractContext, t::Bool) = "bool"

type_string(context::AbstractContext, t::Observable) = type_string(context, t[])

type_string(context::AbstractContext, t::Type) = error("Type $t not supported")

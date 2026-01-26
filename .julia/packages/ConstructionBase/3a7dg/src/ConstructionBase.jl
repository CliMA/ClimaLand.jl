module ConstructionBase

export getproperties
export setproperties
export constructorof
export getfields

# Use markdown files as docstring:
for (name, path) in [
    :ConstructionBase => joinpath(dirname(@__DIR__), "README.md"),
    :constructorof => joinpath(@__DIR__, "constructorof.md"),
    :getfields => joinpath(@__DIR__, "getfields.md"),
    :getproperties => joinpath(@__DIR__, "getproperties.md"),
    :setproperties => joinpath(@__DIR__, "setproperties.md"),
]
    # Don't fail when somehow importing docstrings doesn't work (we
    # don't loose any functionalities for that).  We use explicit
    # `@docs` in `docs/src/index.md` to make sure importing docstrings
    # succeeds in CI.
    try
        include_dependency(path)
        str = read(path, String)
        @eval @doc $str $name
    catch err
        @error "Failed to import docstring for $name" exception=(err, catch_backtrace())
    end
end

constructorof(T::Type) = Base.typename(T).wrapper
constructorof(::Type{<:Tuple}) = tuple
constructorof(::Type{<:NamedTuple{names}}) where names =
    NamedTupleConstructor{names}()

struct NamedTupleConstructor{names} end

@inline function (::NamedTupleConstructor{names})(args...) where names
    NamedTuple{names}(args)
end

################################################################################
#### getfields
################################################################################
getfields(x::Tuple) = x
getfields(x::NamedTuple) = x
getproperties(o::NamedTuple) = o
getproperties(o::Tuple) = o

function check_properties_are_fields(obj)
    # for ntuples of symbols `===` is semantically the same as `==`
    # but triple equals is easier for the compiler to optimize, see #82
    if propertynames(obj) !== fieldnames(typeof(obj))
        error("""
        The `$(nameof(typeof(obj)))` type defines custom properties: it has `propertynames` overloaded.
        Please define `ConstructionBase.setproperties(::$(nameof(typeof(obj))), ::NamedTuple)` to set its properties.
        """)
    end
end

# dispatch on eltype(names) to select Tuple vs NamedTuple
@inline tuple_or_ntuple(names, vals) = tuple_or_ntuple(eltype(names), names, vals)
# if names are empty (object has no properties): return namedtuple, for backwards compat and generally makes more sense than tuple
@inline tuple_or_ntuple(names::Tuple{}, vals::Tuple) = NamedTuple{names}(vals)

# names are consecutive integers: return tuple
@inline function tuple_or_ntuple(::Type{Int}, names, vals)
    @assert Tuple(names) == ntuple(identity, length(names))
    Tuple(vals)
end
# names are symbols: return namedtuple
@inline tuple_or_ntuple(::Type{Symbol}, names, vals::Tuple) = namedtuple(names, vals...)
@inline tuple_or_ntuple(::Type{Symbol}, names, vals) = NamedTuple{Tuple(names)}(vals)
@inline namedtuple(names, vals...) = NamedTuple{Tuple(names)}((vals...,)) # this seemingly unnecessary method encourages union splitting.
# otherwise: throw an error
tuple_or_ntuple(::Type, names, vals) = error("Only Int and Symbol propertynames are supported")

function getproperties(obj)
    fnames = propertynames(obj)
    tuple_or_ntuple(fnames, getproperty.((obj,), fnames))
end
function getfields(obj::T) where {T}
    fnames = fieldnames(T)
    NamedTuple{fnames}(getfield.((obj,), fnames))
end

################################################################################
##### setproperties
################################################################################
function setproperties(obj; kw...)
    setproperties(obj, (;kw...))
end

setproperties(obj             , patch::Tuple      ) = setproperties_object(obj     , patch )
setproperties(obj             , patch::NamedTuple ) = setproperties_object(obj     , patch )
setproperties(obj::Tuple      , patch::Tuple      ) = setproperties_tuple(obj      , patch )
setproperties(obj::Tuple      , patch::NamedTuple ) = setproperties_tuple(obj      , patch )

@generated function check_patch_fields_exist(obj, patch)
    fnames = fieldnames(obj)
    pnames = fieldnames(patch)
    pnames ⊆ fnames ? :(nothing) : :(throw(ArgumentError($("Failed to assign fields $pnames to object with fields $fnames."))))
end

function setproperties(obj::NamedTuple{fields}, patch::NamedTuple{fields}) where {fields}
    patch
end

setproperties_tuple(obj::Tuple, patch::NamedTuple{()}) = obj
@noinline function setproperties_tuple(obj::Tuple, patch::NamedTuple)
    msg = """
    setproperties(obj::Tuple, patch::NamedTuple) only allowed for empty NamedTuple. Got:
    obj  ::Tuple      = $obj
    patch::NamedTuple = $patch
    """
    throw(ArgumentError(msg))
end
function setproperties_tuple(obj::NTuple{N,Any}, patch::NTuple{N,Any}) where {N}
    patch
end
append(x,y) = (x..., y...)
@noinline function throw_setproperties_tuple_error(obj, patch)
    msg = """
    Cannot call `setproperties(obj::Tuple, patch::Tuple)` with `length(obj) < length(patch)`. Got:
    obj = $obj
    patch = $patch
    """
    throw(ArgumentError(msg))
end
function setproperties_tuple(obj::NTuple{N,Any}, patch::NTuple{K,Any}) where {N,K}
    if K > N
        throw_setproperties_tuple_error(obj, patch)
    end
    append(patch, after(obj, Val{K}()))
end
function after(xs::Tuple, ::Val{N}) where {N}
    after(Base.tail(xs), Val{N-1}())
end
function after(x::Tuple, ::Val{0})
    x
end

setproperties_object(obj, patch::Tuple{}) = obj
@noinline function setproperties_object(obj, patch::Tuple)
    msg = """
    setproperties(obj, patch::Tuple) only allowed for empty Tuple. Got:
    obj = $obj
    patch = $patch
    """
    throw(ArgumentError(msg))
end
setproperties_object(obj, patch::NamedTuple{()}) = obj

@generated function setfields_object(obj, patch::NamedTuple)
    args = Expr[]
    pnames = fieldnames(patch)
    for fname in fieldnames(obj)
        source = fname in pnames ? :patch : :obj
        push!(args, :(getproperty($source, $(QuoteNode(fname)))))
    end
    :(constructorof(typeof(obj))($(args...)))
end

function setproperties_object(obj, patch::NamedTuple)
    check_properties_are_fields(obj)
    check_patch_fields_exist(obj, patch)
    setfields_object(obj, patch)
end

include("nonstandard.jl")
include("functions.jl")

end # module

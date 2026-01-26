# We define a custom subtype of `AbstractIrrational` and
# define methods for it that are not generalized yet to `AbstractIrrational`
# in https://github.com/JuliaLang/julia/blob/e536c77f4dc693aafc48af910b4fd86b487e900d/base/irrationals.jl

abstract type IrrationalConstant <: AbstractIrrational end

# TODO: Remove definitions if they become available for `AbstractIrrational` in Base
# Ref https://github.com/JuliaLang/julia/pull/48768
function Base.show(io::IO, ::MIME"text/plain", x::IrrationalConstant)
    if get(io, :compact, false)::Bool
        print(io, x)
    else
        print(io, x, " = ", string(float(x))[1:min(end,15)], "...")
    end
end

Base.:(==)(::T, ::T) where {T<:IrrationalConstant} = true
Base.:<(::T, ::T) where {T<:IrrationalConstant} = false
Base.:<=(::T, ::T) where {T<:IrrationalConstant} = true
Base.hash(x::IrrationalConstant, h::UInt) = 3*objectid(x) - h
Base.round(x::IrrationalConstant, r::RoundingMode) = round(float(x), r)

# https://github.com/JuliaLang/julia/pull/50894
Base.typemin(::Type{T}) where {T<:IrrationalConstant} = T()
Base.typemax(::Type{T}) where {T<:IrrationalConstant} = T()

"""
    @irrational sym [val] def [T]

Define an instance named `sym` of a new singleton type `T` representing an irrational constant as subtype of
`IrrationalConstants.IrrationalConstant <: AbstractIrrational`,
with arbitrary-precision definition in terms of `BigFloat`s given by the expression `def`.

Optionally provide a pre-computed `Float64` value `val` which must equal `Float64(def)`.
It will be computed automatically if omitted.

As default, `T` is set to `sym` with the first character converted to uppercase.

An `AssertionError` is thrown when either `big(def) isa BigFloat` or `Float64(val) == Float64(def)`
returns `false`.

# Examples

```jldoctest; setup = :(import IrrationalConstants)
julia> IrrationalConstants.@irrational twoπ 2*big(π)

julia> twoπ
twoπ = 6.2831853071795...

julia> typeof(twoπ)
Twoπ

julia> IrrationalConstants.@irrational sqrt2 1.4142135623730950488 √big(2)

julia> sqrt2
sqrt2 = 1.4142135623730...

julia> IrrationalConstants.@irrational halfτ 3.14159265358979323846 pi

julia> halfτ
halfτ = 3.1415926535897...

julia> IrrationalConstants.@irrational sqrt3 1.7320508075688772 big(3)
ERROR: AssertionError: big($(Expr(:escape, :sqrt3))) isa BigFloat

julia> IrrationalConstants.@irrational sqrt5 2.2360679775 √big(5)
ERROR: AssertionError: Float64($(Expr(:escape, :sqrt5))) == Float64(big($(Expr(:escape, :sqrt5))))
```
"""
macro irrational(sym::Symbol, val::Float64, def::Union{Symbol,Expr}, T::Symbol=Symbol(uppercasefirst(string(sym))))
    irrational(__module__, sym, val, def, T)
end
macro irrational(sym::Symbol, def::Union{Symbol,Expr}, T::Symbol=Symbol(uppercasefirst(string(sym))))
    irrational(__module__, sym, :(big($(esc(sym)))), def, T)
end
function irrational(mod::Module, sym::Symbol, val::Union{Float64,Expr}, def::Union{Symbol,Expr}, T::Symbol)
    if isdefined(mod, T)
        throw(ArgumentError(LazyString("Type `", T, "` of irrational constant `", sym, "` is already defined in module `", mod, "`.")))
    end
    if sym == T
        throw(ArgumentError(LazyString("The name of the irrational constant (", sym, ") and its type (", T, ") cannot be the same. Please choose a different name for the constant or specify a different type name as the last argument to the macro.")))
    end
    esym = esc(sym)
    qsym = esc(Expr(:quote, sym))
    eT = esc(T)
    bigconvert = if isa(def, Symbol)
        # support older Julia versions prior to https://github.com/JuliaLang/julia/pull/51362
        r = VERSION < v"1.12.0-DEV.78" ? :(Base.MPFR.ROUNDING_MODE[]) : :(Base.Rounding.rounding_raw(BigFloat))
        quote
            function Base.BigFloat(::$eT, r::Base.MPFR.MPFRRoundingMode=$r; precision=precision(BigFloat))
                c = BigFloat(; precision=precision)
                ccall(($(string("mpfr_const_", def)), :libmpfr),
                      Cint, (Ref{BigFloat}, Base.MPFR.MPFRRoundingMode), c, r)
                return c
            end
        end
    else
        quote
            function Base.BigFloat(::$eT; precision=precision(BigFloat))
                setprecision(BigFloat, precision) do
                    $(esc(def))
                end
            end
        end
    end
    quote
        struct $T <: IrrationalConstant end
        const $esym = $eT()
        $bigconvert
        let v = $val, v64 = Float64(v), v32 = Float32(v)
            Base.Float64(::$eT) = v64
            Base.Float32(::$eT) = v32
        end
        Base.show(io::IO, ::$eT) = print(io, $qsym)
        @assert isa(big($esym), BigFloat)
        @assert Float64($esym) == Float64(big($esym))
        @assert Float32($esym) == Float32(big($esym))

        # https://github.com/JuliaLang/julia/pull/55886 removed these effects for generic `AbstractIrrational` subtypes
        @static if isdefined(Base, :_throw_argument_error_irrational_to_rational_bigint)
            function Base.Rational{$BigInt}(::$eT)
                Base._throw_argument_error_irrational_to_rational_bigint()
            end
        end
        @static if isdefined(Base, :_irrational_to_rational)
            Base.@assume_effects :foldable function Base.Rational{T}(x::$eT) where {T<:Integer}
                Base._irrational_to_rational(T, x)
            end
        end
        @static if isdefined(Base, :_irrational_to_float)
            Base.@assume_effects :foldable function (::Type{T})(x::$eT, r::RoundingMode) where {T<:Union{Float32,Float64}}
                Base._irrational_to_float(T, x, r)
            end
        end
        @static if isdefined(Base, :_rationalize_irrational)
            Base.@assume_effects :foldable function Base.rationalize(::Type{T}, x::$eT; tol::Real=0) where {T<:Integer}
                Base._rationalize_irrational(T, x, tol)
            end
        end
        @static if isdefined(Base, :_lessrational)
            Base.@assume_effects :foldable function Base.lessrational(rx::Rational, x::$eT)
                Base._lessrational(rx, x)
            end
        end
    end
end

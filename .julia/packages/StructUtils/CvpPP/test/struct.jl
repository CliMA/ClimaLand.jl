struct A
    a::Int
    b::Int
    c::Int
    d::Int
end

struct AA
    a::Int
    b::Int
    c::Int
    d::Int
    e::Int
end

StructUtils.fielddefaults(::StructUtils.StructStyle, ::Type{AA}) = (e=5,)

mutable struct B
    a::Int
    b::Int
    c::Int
    d::Int
    B() = new()
end

StructUtils.noarg(::StructUtils.StructStyle, ::Type{B}) = true
Base.:(==)(b1::B, b2::B) = b1.a == b2.a && b1.b == b2.b && b1.c == b2.c && b1.d == b2.d

Base.@kwdef struct BB
    a::Int = 1
    b::Int = 2
    c::Int = 3
    d::Int = 4
end

StructUtils.kwarg(::StructUtils.StructStyle, ::Type{BB}) = true

struct C
end

struct D
    a::Int
    b::Float64
    c::String
end

struct LotsOfFields
    x1::String
    x2::String
    x3::String
    x4::String
    x5::String
    x6::String
    x7::String
    x8::String
    x9::String
    x10::String
    x11::String
    x12::String
    x13::String
    x14::String
    x15::String
    x16::String
    x17::String
    x18::String
    x19::String
    x20::String
    x21::String
    x22::String
    x23::String
    x24::String
    x25::String
    x26::String
    x27::String
    x28::String
    x29::String
    x30::String
    x31::String
    x32::String
    x33::String
    x34::String
    x35::String
end

struct Wrapper
    x::NamedTuple{(:a, :b), Tuple{Int, String}}
end

mutable struct UndefGuy
    id::Int
    name::String
    UndefGuy() = new()
end

StructUtils.noarg(::StructUtils.StructStyle, ::Type{UndefGuy}) = true

struct E
    id::Int
    a::A
end

Base.@kwdef struct F
    id::Int
    rate::Float64
    name::String
end

StructUtils.kwarg(::StructUtils.StructStyle, ::Type{F}) = true

Base.@kwdef struct G
    id::Int
    rate::Float64
    name::String
    f::F
end

StructUtils.kwarg(::StructUtils.StructStyle, ::Type{G}) = true

struct H
    id::Int
    name::String
    properties::Dict{String, Any}
    addresses::Vector{String}
end

@enum Fruit apple banana

struct I
    id::Int
    name::String
    fruit::Fruit
end

abstract type Vehicle end

struct Car <: Vehicle
    make::String
    model::String
    seatingCapacity::Int
    topSpeed::Float64
end

struct Truck <: Vehicle
    make::String
    model::String
    payloadCapacity::Float64
end

struct J
    id::Union{Int, Nothing}
    name::Union{String, Nothing}
    rate::Union{Int, Float64}
end

struct K
    id::Int
    value::Union{Float64, Missing}
end

Base.@kwdef struct System
    duration::Real = 0 # mandatory
    cwd::Union{Nothing, String} = nothing
    environment::Union{Nothing, Dict} = nothing
    batch::Union{Nothing, Dict} = nothing
    shell::Union{Nothing, Dict} = nothing
end

StructUtils.kwarg(::StructUtils.StructStyle, ::Type{System}) = true

struct L
    id::Int
    first_name::String
    rate::Float64
end

struct ThreeDates
    date::Date
    datetime::DateTime
    time::Time
end

struct M
    id::Int
    value::Union{Nothing,K}
end

struct Recurs
    id::Int
    value::Union{Nothing,Recurs}
end

struct N
    id::Int
    uuid::UUID
end

@tags struct O
    id::Int
    name::Union{I,L,Missing,Nothing} &(choosetype=x->isnothing(x) ? Nothing : ismissing(x) ? Missing : haskey(x, :fruit) ? I : L,)
end

@noarg mutable struct P
    id::Int
    @atomic(name::String) = "Jim"
end

struct Point
    x::Int
    y::Int
end

# Test structs for @nonstruct macro
@nonstruct struct NonStructUnit
    value::String
end

@nonstruct struct NonStructComplex
    id::Int
    data::String
end

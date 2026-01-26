using ConstructionBase
using Test
using LinearAlgebra

struct Empty end
struct AB{A,B}
    a::A
    b::B
end

@testset "constructorof" begin
    @test constructorof(Empty)() === Empty()
    @inferred constructorof(AB{Int, Int})
    @test constructorof(AB{Int, Int})(1, 2) === AB(1,2)
    @test constructorof(AB{Int, Int})(1.0, 2) === AB(1.0,2)
    @test constructorof(typeof((a=1, b=2)))(1.0, 2) === (a=1.0, b=2)
    @test constructorof(NamedTuple{(:a, :b)})(1.0, 2) === (a=1.0, b=2)
    @test constructorof(Tuple)(1.0, 2) === (1.0, 2)
    @test constructorof(Tuple{Nothing, Missing})(1.0, 2) === (1.0, 2)
end

@testset "getfields" begin
    @test getfields(()) === ()
    @test keys(getfields([])) == fieldnames(typeof([]))
    @test getfields(Empty()) === NamedTuple()
    @test getfields(NamedTuple()) === NamedTuple()
    @test getfields((10,20,30)) === (10,20,30)
    @test getfields((a=10,b=20f0,c=true)) === (a=10,b=20f0,c=true)
    @test getfields(AB(1, 10)) === (a=1, b=10)
    adder(a) = x -> x + a
    @test getfields(adder(1)) === (a=1,)
end

struct DontTouchProperties
    a
    b
end
Base.propertynames(::DontTouchProperties) = error()
Base.getproperty(::DontTouchProperties, ::Symbol) = error()
ConstructionBase.getproperties(::DontTouchProperties) = error()
@testset "getfields does not depend on properties" begin
    @test getfields(DontTouchProperties(1,2)) === (a=1, b=2)
    @test constructorof(DontTouchProperties) === DontTouchProperties
end

@testset "getproperties" begin
    o = AB(1, 2)
    @test getproperties(o) === (a=1, b=2)
    @inferred getproperties(o)
    @test getproperties(Empty()) === NamedTuple()

    t = (1,2,3.0)
    @inferred getproperties(t)
    @test t === getproperties(t)
    @test () === getproperties(())
end

@testset "setproperties" begin

    @test setproperties(NamedTuple(), NamedTuple()) === NamedTuple()
    @test setproperties((), NamedTuple()) === ()
    @test setproperties(NamedTuple(), ()) === NamedTuple()
    @test setproperties((), ()) === ()
    @test setproperties(1, ()) === 1
    @test setproperties(1, NamedTuple()) === 1

    @test setproperties((1,), ()) === (1,)
    @test setproperties((1,), NamedTuple()) === (1,)
    @test setproperties((a=1,), ()) === (a=1,)
    @test setproperties((a=1,), NamedTuple()) === (a=1,)
    @test setproperties(AB(1,2), ()) === AB(1,2)
    @test_throws ArgumentError setproperties(AB(1,2), (10,))
    @test setproperties(AB(1,2), NamedTuple()) === AB(1,2)

    @test setproperties(AB(1,2), (a=2, b=3))   === AB(2,3)
    @test setproperties(AB(1,2), (a=2, b=3.0)) === AB(2,3.0)
    @test setproperties(AB(1,2), a=2, b=3.0) === AB(2,3.0)

    res = @test_throws ArgumentError setproperties(AB(1,2), (a=2, this_field_does_not_exist=3.0))
    msg = sprint(showerror, res.value)
    @test occursin("this_field_does_not_exist", msg)

    res = @test_throws ArgumentError setproperties(AB(1,2), a=2, this_field_does_not_exist=3.0)
    msg = sprint(showerror, res.value)
    @test occursin("this_field_does_not_exist", msg)

    @test setproperties(42, NamedTuple()) === 42
    @test setproperties(42) === 42

    @test setproperties(Empty(), NamedTuple()) === Empty()
    @test setproperties(Empty()) === Empty()

    @test setproperties((a=1, b=2), (a=1.0,)) === (a=1.0, b=2)
    @test setproperties((a=1, b=2), a=1.0) === (a=1.0, b=2)

    @inferred setproperties(AB(1,2), a=2, b=3.0)
    @inferred setproperties(Empty(), NamedTuple())
    @inferred setproperties((a=1, b=2), a=1.0)
    @inferred setproperties((a=1, b=2), (a=1.0,))

    @test setproperties((1,), ()) === (1,)
    @test setproperties((1,), (10,)) === (10,)
    @test_throws ArgumentError setproperties((1,), (10,20)) === (10,)

    @inferred setproperties((1,2,3), (1,2,3))
    @test setproperties((1,2,3), ()) === (1,2,3)
    @test setproperties((1,2,3), (10.0,)) === (10.0,2,3)
    @test setproperties((1,2,3), (10.0,20)) === (10.0,20,3)
    @test setproperties((1,2,3), (10.0,20,30)) === (10.0,20,30)
    @test_throws ArgumentError setproperties((1,2,3), (10.0,20,30,40))

    @test_throws ArgumentError setproperties((a=1,b=2), (10,20))
    @test_throws ArgumentError setproperties((), (10,))
    @test_throws ArgumentError setproperties((1,2), (a=10,b=20))
end

struct CustomSetproperties
    _a::Int
end
function ConstructionBase.setproperties(o::CustomSetproperties, patch::NamedTuple)
    if isempty(patch)
        o
    elseif propertynames(patch) == (:a,)
        CustomSetproperties(patch.a)
    else
        error()
    end
end

@testset "custom setproperties unambiguous on empty" begin
    o = CustomSetproperties(1)
    @test o === setproperties(o)
    @test o === setproperties(o, NamedTuple())
    @test CustomSetproperties(2) === setproperties(o, a=2)
    @test CustomSetproperties(2) === setproperties(o, (a=2,))
end

@testset "constructors for non-standadard Base and LinearAlgebra etc objects" begin
    A1 = zeros(5, 6)
    A2 = ones(Float32, 5, 6)

    @testset "SubArray" begin
        subarray = view(A1, 1:2, 3:4)
        @test constructorof(typeof(subarray))(getfields(subarray)...) === subarray
        @test all(constructorof(typeof(subarray))(A2, (Base.OneTo(2), 3:4), 0, 0) .== Float32[1 1; 1 1])
        @inferred constructorof(typeof(subarray))(getfields(subarray)...)
        @inferred constructorof(typeof(subarray))(A2, (Base.OneTo(2), 3:4), 0, 0)
    end

    @testset "ReinterpretArray" begin
        ra1 = reinterpret(Float16, A1)
        @test constructorof(typeof(ra1))(A1) === ra1
        @test constructorof(typeof(ra1))(getfields(ra1)...) === ra1
        ra2 = constructorof(typeof(ra1))(A2)
        @test size(ra2) == (10, 6)
        @test eltype(ra2) == Float16
        @inferred constructorof(typeof(ra1))(getfields(ra1)...)
        @inferred constructorof(typeof(ra1))(A2)
    end

    @testset "PermutedDimsArray" begin
        pda1 = PermutedDimsArray(A1, (2, 1))
        @test constructorof(typeof(pda1))(A1) === pda1
        @test constructorof(typeof(pda1))(getfields(pda1)...) === pda1
        @test eltype(constructorof(typeof(pda1))(A2)) == Float32
        @inferred constructorof(typeof(pda1))(getfields(pda1)...)
        @inferred constructorof(typeof(pda1))(A2)
    end

    @testset "Tridiagonal" begin
        d = randn(12)
        dl = randn(11)
        du = randn(11)
        tda = Tridiagonal(dl, d, du)
        @test isdefined(tda, :du2) == false
        @test constructorof(typeof(tda))(dl, d, du) === tda
        @test constructorof(typeof(tda))(getfields(tda)...) === tda
        # lu factorization defines du2
        tda_lu = lu!(tda).factors
        @test isdefined(tda_lu, :du2) == true
        @test constructorof(typeof(tda_lu))(getfields(tda_lu)...) === tda_lu
        @test constructorof(typeof(tda_lu))(getfields(tda)...) !== tda_lu
        @test constructorof(typeof(tda_lu))(getfields(tda)...) === tda
        @inferred constructorof(typeof(tda))(getfields(tda)...)
        @inferred constructorof(typeof(tda))(getfields(tda_lu)...)
    end

    @testset "LinRange" begin
        lr1 = LinRange(1, 7, 10)
        lr2 = LinRange(1.0f0, 7.0f0, 10)
        @test constructorof(typeof(lr1))(1, 7, 10, nothing) === lr1
        @test constructorof(typeof(lr1))(getfields(lr2)...) === lr2
        @inferred constructorof(typeof(lr1))(getfields(lr1)...)
        @inferred constructorof(typeof(lr1))(getfields(lr2)...)
    end

    @testset "Cholesky" begin
        x = randn(3, 3)
        X = x * x'
        @testset "uplo=$uplo" for uplo in ['L', 'U']
            C = Cholesky(X, uplo, 0)

            # Empty patch.
            C_new = ConstructionBase.setproperties(C, NamedTuple())
            @test typeof(C_new) === typeof(C)
            for f in propertynames(C)
                @test getproperty(C_new, f) == getproperty(C, f)
            end

            # Update `L`.
            C_new = ConstructionBase.setproperties(C, (L=2 * C.L,))
            @test typeof(C_new) === typeof(C)
            for f in propertynames(C)
                @test getproperty(C_new, f) == 2 * getproperty(C, f)
            end

            # Update `U`.
            C_new = ConstructionBase.setproperties(C, (U=2 * C.U,))
            @test typeof(C_new) === typeof(C)
            for f in propertynames(C)
                @test getproperty(C_new, f) == 2 * getproperty(C, f)
            end

            # Update `UL`
            C_new = ConstructionBase.setproperties(C, (UL=2 * C.UL,))
            @test typeof(C_new) === typeof(C)
            for f in propertynames(C)
                @test getproperty(C_new, f) == 2 * getproperty(C, f)
            end

            # We can only set the properties with `LowerTriangular` or `UpperTriangular` matrices.
            @test_throws ArgumentError ConstructionBase.setproperties(C, (L=parent(C.L),))
            @test_throws ArgumentError ConstructionBase.setproperties(C, (U=parent(C.U),))
            # Can only set one at the time.
            @test_throws ArgumentError ConstructionBase.setproperties(C, (L=C.L, U=C.U,))
            @test_throws ArgumentError ConstructionBase.setproperties(C, (UL=C.UL, U=C.U,))
            @test_throws ArgumentError ConstructionBase.setproperties(C, (UL=C.UL, L=C.L,))
            # And make sure any other patch will fail.
            @test_throws ArgumentError ConstructionBase.setproperties(C, (asdf=C.UL,))
        end
    end

    @testset "Expr" begin
        e = :(a + b)
        @test e == @inferred constructorof(typeof(e))(getfields(e)...)
    end
end

@testset "default function constructors" begin
    function multiplyer(a, b)
        x -> x * a * b
    end

    mult11 = multiplyer(1, 1)
    @test mult11(1) === 1
    mult23 = @inferred constructorof(typeof(mult11))(2.0, 3.0)
    @inferred mult23(1)
    @test mult23(1) === 6.0
    multbc = @inferred constructorof(typeof(mult23))("b", "c")
    @inferred multbc("a")
    @test multbc("a") == "abc"

    @test (@inferred constructorof(typeof(sin))()) === sin
    f = x -> x^2
    @test (@inferred constructorof(typeof(f))()) === f
end

struct Adder{V} <: Function
    value::V
end
(o::Adder)(x) = o.value + x

struct Adder2{V} <: Function
    value::V
    int::Int
end
(o::Adder2)(x) = o.value + o.int + x

struct AddTuple{T} <: Function
    tuple::Tuple{T,T,T}
end
(o::AddTuple)(x) = sum(o.tuple) + x

# A function with an inner constructor with checks
struct Rotation{M} <: Function
    matrix::M
    function Rotation(m)
        @assert isapprox(det(m), 1)
        @assert isapprox(m*m', I)
        new{typeof(m)}(m)
    end
end

@testset "Custom function object constructors still work" begin
    add1 = Adder(1)
    @test add1(1) === 2
    add2 = @inferred ConstructionBase.constructorof(typeof(add1))(2.0)
    @inferred add2(1)
    @test add2(1) == 3.0
    add12 = Adder2(1, 2)
    @test @inferred add12(3) ==  6
    add22 = @inferred ConstructionBase.constructorof(typeof(add12))(2.0, 2)
    @test @inferred add22(3) ==  7.0

    addtuple123 = AddTuple((1, 2, 3))
    @test addtuple123(1) === 7
    addtuple234 = @inferred ConstructionBase.constructorof(typeof(addtuple123))((2.0, 3.0, 4.0))
    @inferred addtuple234(1)
    @test addtuple234(1) === 10.0

    @testset "inner constructor without type parameters is still called" begin
        @test_throws AssertionError constructorof(Rotation{Matrix{Float64}})(zeros(3,3))
    end
end

# example of a struct with different fields and properties
struct FieldProps{NT <: NamedTuple{(:a, :b)}}
    components::NT
end

Base.propertynames(::FieldProps) = (:a, :b)
Base.getproperty(obj::FieldProps, name::Symbol) = getproperty(getfield(obj, :components), name)

@testset "use properties, not fields" begin
    x = FieldProps((a=1, b=:b))
    @test constructorof(typeof(x)) === FieldProps
    @test getfields(x) === (components=(a=1, b=:b),)
    res = @test_throws ErrorException setproperties(x, c=0)
    msg = sprint(showerror, res.value)
    @test occursin("overload", msg)
    @test occursin("setproperties", msg)
    @test occursin("FieldProps", msg)
    @test_throws ErrorException setproperties(x, components=(a=1,b=:b))
    msg = sprint(showerror, res.value)
    @test occursin("overload", msg)
    @test occursin("setproperties", msg)
    @test occursin("FieldProps", msg)
    @test_throws ErrorException setproperties(x, a="aaa")
    msg = sprint(showerror, res.value)
    @test occursin("overload", msg)
    @test occursin("setproperties", msg)
    @test occursin("FieldProps", msg)
 # == FieldProps((a="aaa", b=:b)
    @test getproperties(x) == (a=1, b=:b)
end


struct SProp
    names
end
Base.propertynames(s::SProp) = getfield(s, :names)
Base.getproperty(s::SProp, prop::Symbol) = "ps$prop"
Base.getproperty(s::SProp, prop::Int) = "pi$prop"
Base.getproperty(s::SProp, prop::String) = "pstr$prop"

@testset "properties can be numbered" begin
    @test getproperties(SProp((:a, :b))) === (a="psa", b="psb")
    @test getproperties(SProp((1, 2))) === ("pi1", "pi2")
    # what should it return?
    @test_broken getproperties(SProp(("a", "b")))

    @test_throws ErrorException getproperties(SProp((1, :a)))
end

@testset "propertynames can be a vector" begin
    @test getproperties(SProp([:a, :b])) === (a="psa", b="psb")
    @test getproperties(SProp(Symbol[])) === (;)
    @test getproperties(SProp([1, 2])) === ("pi1", "pi2")
    @test getproperties(SProp(Int[])) === ()
end

function funny_numbers(::Type{Tuple}, n)::Tuple
    types = [
        Int128, Int16, Int32, Int64, Int8,
        UInt128, UInt16, UInt32, UInt64, UInt8,
        Float16, Float32, Float64,
    ]
    Tuple([T(true) for T in rand(types, n)])
end

function funny_numbers(::Type{NamedTuple}, n)::NamedTuple
    t = funny_numbers(Tuple,n)
    pairs = map(1:n) do i
        Symbol("a$i") => t[i]
    end
    (;pairs...)
end

abstract type S end
Sn_from_n = Dict{Int,Type}()
for n in [0,1,10,20,40]
    Sn = Symbol("S$n")
    types = [Symbol("A$i") for i in 1:n]
    fields = [Symbol("a$i") for i in 1:n]
    typed_fields = [:($ai::$Ai) for (ai,Ai) in zip(fields, types)]
    @eval struct $(Sn){$(types...)} <: S
        $(typed_fields...)
    end
    @eval Sn_from_n[$n] = $Sn
end
function funny_numbers(::Type{S}, n)::S
    fields = funny_numbers(Tuple, n)
    Sn_from_n[n](fields...)
end

reconstruct(obj, content) = constructorof(typeof(obj))(content...)

function write_output_to_ref!(f, out_ref::Ref, arg_ref::Ref)
    arg = arg_ref[]
    out_ref[] = f(arg)
    out_ref
end
function write_output_to_ref!(f, out_ref::Ref, arg_ref1::Ref, arg_ref2::Ref)
    arg1 = arg_ref1[]
    arg2 = arg_ref2[]
    out_ref[] = f(arg1,arg2)
    out_ref
end
function hot_loop_allocs(f::F, args...) where {F}
    # we want to test that f(args...) does not allocate
    # when used in hot loops
    # however a naive @allocated f(args...)
    # will not be representative of what happens in an inner loop
    # Instead it will sometimes box inputs/outputs
    # and report too many allocations
    # so we use Refs to minimize inputs and outputs
    out_ref = Ref(f(args...))
    arg_refs = map(Ref, args)
    write_output_to_ref!(f, out_ref, arg_refs...)
    out_ref = typeof(out_ref)() # erase out_ref so we can assert work was done later
    # Avoid splatting args... which also results in undesired allocs
    allocs = if length(arg_refs) == 1
        r1, = arg_refs
        @allocated write_output_to_ref!(f, out_ref, r1)
    elseif length(arg_refs) == 2
        r1,r2 = arg_refs
        @allocated write_output_to_ref!(f, out_ref, r1, r2)
    else
        error("TODO too many args")
    end
    @assert out_ref[] == f(args...)
    return allocs
end

@testset "no allocs $T" for T in [Tuple, NamedTuple, S]
    @testset "n = $n" for n in [0,1,10,20]
        obj = funny_numbers(T, n)
        new_content = funny_numbers(Tuple, n)
        @test 0 == hot_loop_allocs(constructorof, typeof(obj))
        @test 0 == hot_loop_allocs(reconstruct, obj, new_content)
        @test 0 == hot_loop_allocs(getproperties, obj)
        @test 0 == hot_loop_allocs(getfields, obj)
        patch_sizes = [0,1,n÷3,n÷2,n]
        patch_sizes = min.(patch_sizes, n)
        patch_sizes = unique(patch_sizes)
        for k in patch_sizes
            patch = if T === Tuple
                funny_numbers(Tuple, k)
            else
                funny_numbers(NamedTuple, k)
            end
            @test 0 == hot_loop_allocs(setproperties, obj, patch)
        end
    end
end

struct S2
    a::Union{Nothing, Int}
    b::Union{UInt32, Int32}
end

@testset "no allocs S2" begin
    obj = S2(3, UInt32(5))
    @test 0 == hot_loop_allocs(constructorof, typeof(obj))
    @test 0 == hot_loop_allocs(setproperties, obj, (; a = nothing, b = Int32(6)))
end

@testset "inference" begin
    @testset "Tuple n=$n" for n in [0,1,2,3,4,5,10,20,30,40]
        t = funny_numbers(Tuple,n)
        @test length(t) == n
        @test getproperties(t) === t
        @inferred getproperties(t)
        @test getfields(t) === t
        @inferred getfields(t)
        @inferred constructorof(typeof(t))
        content = funny_numbers(Tuple,n)
        @inferred reconstruct(t, content)
        for k in 0:n
            t2 = funny_numbers(Tuple,k)
            @test setproperties(t, t2)[1:k] === t2
            @test setproperties(t, t2) isa Tuple
            @test length(setproperties(t, t2)) == n
            @test setproperties(t, t2)[k+1:n] === t[k+1:n]
            @inferred setproperties(t, t2)
        end
    end
    @inferred getproperties(funny_numbers(Tuple,100))
    @inferred setproperties(funny_numbers(Tuple,100), funny_numbers(Tuple,90))

    @testset "NamedTuple n=$n" for n in [0,1,2,3,4,5,10,20,30,40]
        nt = funny_numbers(NamedTuple, n)
        @test nt isa NamedTuple
        @test length(nt) == n
        @test getproperties(nt) === nt
        @inferred getproperties(nt)
        @test getfields(nt) === nt
        @inferred getfields(nt)

        @inferred constructorof(typeof(nt))
        content = funny_numbers(NamedTuple,n)
        @inferred reconstruct(nt, content)
        #no_allocs_test(nt, content)
        for k in 0:n
            nt2 = funny_numbers(NamedTuple, k)
            @inferred setproperties(nt, nt2)
            @test Tuple(setproperties(nt, nt2))[1:k] === Tuple(nt2)
            @test setproperties(nt, nt2) isa NamedTuple
            @test length(setproperties(nt, nt2)) == n
            @test Tuple(setproperties(nt, nt2))[k+1:n] === Tuple(nt)[k+1:n]
        end
    end
    @inferred getproperties(funny_numbers(NamedTuple, 100))
    @inferred setproperties(funny_numbers(NamedTuple, 100), funny_numbers(NamedTuple, 90))

    @inferred setproperties(funny_numbers(S,0), funny_numbers(NamedTuple, 0))
    @inferred setproperties(funny_numbers(S,1), funny_numbers(NamedTuple, 1))
    @inferred setproperties(funny_numbers(S,20), funny_numbers(NamedTuple, 18))
    @inferred setproperties(funny_numbers(S,40), funny_numbers(NamedTuple, 38))
    @inferred constructorof(S0)
    @inferred constructorof(S1)
    @inferred constructorof(S20)
    @inferred constructorof(S40)
    @inferred reconstruct(funny_numbers(S,0) , funny_numbers(Tuple,0))
    @inferred reconstruct(funny_numbers(S,1) , funny_numbers(Tuple,1))
    @inferred reconstruct(funny_numbers(S,20), funny_numbers(Tuple,20))
    @inferred reconstruct(funny_numbers(S,40), funny_numbers(Tuple,40))

    @inferred getfields(funny_numbers(S,0))
    @inferred getfields(funny_numbers(S,1))
    @inferred getfields(funny_numbers(S,20))
    @inferred getfields(funny_numbers(S,40))
    @inferred getproperties(funny_numbers(S,0))
    @inferred getproperties(funny_numbers(S,1))
    @inferred getproperties(funny_numbers(S,20))
    @inferred getproperties(funny_numbers(S,40))
end


using StaticArrays, IntervalSets

@testset "staticarrays" begin
    sa = @SVector [2, 4, 6, 8]
    sa2 = ConstructionBase.constructorof(typeof(sa))((3.0, 5.0, 7.0, 9.0))
    @test sa2 === @SVector [3.0, 5.0, 7.0, 9.0]

    ma = @MMatrix [2.0 4.0; 6.0 8.0]
    ma2 = @inferred ConstructionBase.constructorof(typeof(ma))((1, 2, 3, 4))
    @test ma2 isa MArray{Tuple{2,2},Int,2,4}
    @test all(ma2 .=== @MMatrix [1 3; 2 4])

    sz = SizedArray{Tuple{2,2}}([1 2;3 4])
    sz2 = @inferred ConstructionBase.constructorof(typeof(sz))([:a :b; :c :d])
    @test sz2 == SizedArray{Tuple{2,2}}([:a :b; :c :d])
    @test typeof(sz2) <: SizedArray{Tuple{2,2},Symbol,2,2}

    for T in (SVector, MVector)
        @test @inferred(ConstructionBase.constructorof(T)((1, 2, 3)))::T == T((1, 2, 3))
        @test @inferred(ConstructionBase.constructorof(T{3})((1, 2, 3)))::T == T((1, 2, 3))
        @test @inferred(ConstructionBase.constructorof(T{3})((1, 2)))::T == T((1, 2))
        @test @inferred(ConstructionBase.constructorof(T{3, Symbol})((1, 2, 3)))::T == T((1, 2, 3))
        @test @inferred(ConstructionBase.constructorof(T{3, Symbol})((1, 2)))::T == T((1, 2))
        @test @inferred(ConstructionBase.constructorof(T{3, X} where {X})((1, 2, 3)))::T == T((1, 2, 3))
        @test @inferred(ConstructionBase.constructorof(T{3, X} where {X})((1, 2)))::T == T((1, 2))
        @test @inferred(ConstructionBase.constructorof(T{X, Symbol} where {X})((1, 2, 3)))::T == T((1, 2, 3))
    end

    sv = SVector(1, 2)
    @test SVector(3.0, 2.0) === @inferred setproperties(sv, x = 3.0)
    @test SVector(3.0, 5.0) === @inferred setproperties(sv, x = 3.0, y = 5.0)
    @test SVector(-1.0, -2.0) === @inferred setproperties(sv, data = (-1.0, -2))
    @test_throws "does not have properties z" setproperties(sv, z = 3.0)
    @test_throws "does not have properties z" setproperties(SVector(1, 2, 3, 4, 5), z = 3.0)
end

@testset "intervalsets" begin
    @test constructorof(typeof(1..2))(0.5, 1.5) === 0.5..1.5
    @test constructorof(typeof(OpenInterval(1, 2)))(0.5, 1.5) === OpenInterval(0.5, 1.5)
    @test setproperties(1..2, left=0.0) === 0.0..2.0
    @test setproperties(OpenInterval(1.0, 2.0), left=1, right=5) === OpenInterval(1, 5)
end

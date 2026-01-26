using Static, Aqua
using Static: Zero
using Test

@testset "Aqua" begin
    Aqua.find_persistent_tasks_deps(Static)
    Aqua.test_ambiguities(Static, recursive = false)
    Aqua.test_deps_compat(Static)
    Aqua.test_piracies(Static, broken = true)
    Aqua.test_project_extras(Static)
    Aqua.test_stale_deps(Static)
    Aqua.test_unbound_args(Static)
    Aqua.test_undefined_exports(Static)
end

@testset "StaticSymbol" begin
    x = StaticSymbol(:x)
    y = StaticSymbol("y")
    z = StaticSymbol(1)
    @test y === StaticSymbol(:y)
    @test z === StaticSymbol(Symbol(1))
    @test z == StaticSymbol(Symbol(1))
    @test @inferred(StaticSymbol(x)) === x
    @test @inferred(StaticSymbol(x, y)) === StaticSymbol(:x, :y)
    @test @inferred(StaticSymbol(x, y, z)) === static(:xy1)
    @test @inferred(static(nothing)) === nothing
    nt = (x = :x, y = "y", z = 1)
    @test @inferred(getindex(nt, x)) == :x
    @test @inferred(getindex(nt, (x, y))) == (x = :x, y = "y")
    @test_throws ErrorException static([])
end

@testset "StaticInt" begin
    @test static(UInt(8)) === StaticInt(UInt(8)) === StaticInt{8}()
    @test iszero(StaticInt(0))
    @test !iszero(StaticInt(1))
    @test !isone(StaticInt(0))
    @test isone(StaticInt(1))
    @test @inferred(one(StaticInt(1))) === StaticInt(1)
    @test @inferred(zero(StaticInt(1))) === StaticInt(0)
    @test @inferred(one(StaticInt)) === StaticInt(1)
    @test @inferred(zero(StaticInt)) === StaticInt(0) === StaticInt(StaticInt(Val(0)))
    @test eltype(one(StaticInt)) <: Int
    @test @inferred(isodd(one(StaticInt)))
    @test @inferred(iseven(zero(StaticInt)))

    x = StaticInt(1)
    y = static(3)
    @test isinteger(y)
    @test @inferred(minmax(x, y)) == @inferred(minmax(y, x)) == minmax(1, 3)
    @test @inferred(Bool(x)) isa Bool
    @test @inferred(BigInt(x)) isa BigInt
    @test @inferred(Integer(x)) === x
    @test @inferred(%(x, Integer)) === 1
    # test for ambiguities and correctness
    for i in Any[StaticInt(0), StaticInt(1), StaticInt(2), 3]
        for j in Any[StaticInt(0), StaticInt(1), StaticInt(2), 3]
            i === j === 3 && continue
            for f in [+, -, *, ÷, %, <<, >>, >>>, &, |, ⊻, ==, ≤, ≥, min, max]
                (iszero(j) && ((f === ÷) || (f === %))) && continue # integer division error
                @test convert(Int, @inferred(f(i, j))) ==
                      f(convert(Int, i), convert(Int, j))
            end
        end
        i == 3 && break
        for f in [+, -, *, /, ÷, %, ==, ≤, ≥]
            w = f(convert(Int, i), 1.4)
            x = f(1.4, convert(Int, i))
            @test convert(typeof(w), @inferred(f(i, 1.4))) === w
            @test convert(typeof(x), @inferred(f(1.4, i))) === x # if f is division and i === StaticInt(0), returns `NaN`; hence use of ==== in check.
            (((f === ÷) || (f === %)) && (i === StaticInt(0))) && continue
            y = f(convert(Int, i), 2 // 7)
            z = f(2 // 7, convert(Int, i))
            @test convert(typeof(y), @inferred(f(i, 2 // 7))) === y
            @test convert(typeof(z), @inferred(f(2 // 7, i))) === z
        end
    end
    @test @inferred(*(Zero(), 3)) === @inferred(*(3, Zero())) === *(Zero(), Zero())

    @test float(StaticInt(8)) === static(8.0)
    @test @inferred(cld(40, static(4))) === cld(40, 4)
    @test @inferred(fld(40, static(4))) === fld(40, 4)
    @test @inferred(cld(static(40), static(4))) === static(cld(40, 4))
    @test @inferred(fld(static(40), static(4))) === static(fld(40, 4))

    # test specific promote rules to ensure we don't cause ambiguities
    SI = StaticInt{1}
    IR = typeof(1 // 1)
    PI = typeof(pi)
    @test @inferred(convert(SI, SI())) === SI()
    @test @inferred(promote_type(SI, PI)) <: promote_type(Int, PI)
    @test @inferred(promote_type(SI, IR)) <: promote_type(Int, IR)
    @test @inferred(promote_rule(SI, SI)) <: Int
    @test @inferred(promote_type(Base.TwicePrecision{Int}, StaticInt{1})) <:
          Base.TwicePrecision{Int}

    @test static(Int8(-18)) === static(-18)
    @test static(0xef) === static(239)
    @test static(Int16(-18)) === static(-18)
    @test static(0xffef) === static(65519)
    if sizeof(Int) == 8
        @test static(Int32(-18)) === static(-18)
        @test static(0xffffffef) === static(4294967279)
    end
    @test @inferred(getindex([1], static(1))) == 1
    @test @inferred(Base.checkindex(Bool, 1:10, static(1)))
    @test widen(static(1)) isa Int128
    @test isless(static(1), static(2))

    v = rand(3)
    @test real(static(3)) === static(3)
    @test imag(static(3)) === static(0)

    @test mod(static(3), static(3)) === static(0)
    @test mod(static(3), 3) == 0
    @test mod(3, static(3)) == 0
    @test mod(static(3), static(2)) === static(1)
    @test mod(static(3), 2) == 1
    @test mod(3, static(2)) == 1
end

@testset "StaticBool" begin
    t = static(static(true))
    f = StaticBool(static(false))

    @test @inferred(dynamic(t)) === @inferred(known(t)) === true
    @test StaticBool{true}() === t
    @test StaticBool{false}() === f

    @test @inferred(StaticInt(t)) === StaticInt(1)
    @test @inferred(StaticInt(f)) === StaticInt(0)

    @test isinteger(static(true))

    @test @inferred(~t) === f
    @test @inferred(~f) === t
    @test @inferred(!t) === f
    @test @inferred(!f) === t
    @test @inferred(+t) === StaticInt(1)
    @test @inferred(+f) === StaticInt(0)
    @test @inferred(-t) === StaticInt(-1)
    @test @inferred(-f) === StaticInt(0)
    @test @inferred(sign(t)) === t
    @test @inferred(abs(t)) === t
    @test @inferred(abs2(t)) === t
    @test @inferred(iszero(zero(True)))
    @test @inferred(isone(one(True)))
    @test !@inferred(iszero(t))
    @test @inferred(isone(t))
    @test @inferred(iszero(f))
    @test !@inferred(isone(f))

    @test @inferred(xor(true, f))
    @test @inferred(xor(f, true))
    @test @inferred(xor(f, f)) === f
    @test @inferred(xor(f, t)) === t
    @test @inferred(xor(t, f)) === t
    @test @inferred(xor(t, t)) === f

    @test @inferred(|(true, f))
    @test @inferred(|(true, t)) === t
    @test @inferred(|(f, true))
    @test @inferred(|(t, true)) === t
    @test @inferred(|(f, f)) === f
    @test @inferred(|(f, t)) === t
    @test @inferred(|(t, f)) === t
    @test @inferred(|(t, t)) === t

    @test @inferred(Base.:(&)(true, f)) === f
    @test @inferred(Base.:(&)(true, t))
    @test @inferred(Base.:(&)(f, true)) === f
    @test @inferred(Base.:(&)(t, true))
    @test @inferred(Base.:(&)(f, f)) === f
    @test @inferred(Base.:(&)(f, t)) === f
    @test @inferred(Base.:(&)(t, f)) === f
    @test @inferred(Base.:(&)(t, t)) === t

    @test @inferred(<(f, f)) === false
    @test @inferred(<(f, t)) === true
    @test @inferred(<(t, f)) === false
    @test @inferred(<(t, t)) === false

    @test @inferred(<=(f, f)) === true
    @test @inferred(<=(f, t)) === true
    @test @inferred(<=(t, f)) === false
    @test @inferred(<=(t, t)) === true

    @test @inferred(==(f, f)) === true
    @test @inferred(==(f, t)) === false
    @test @inferred(==(t, f)) === false
    @test @inferred(==(t, t)) === true

    @test @inferred(*(f, t)) === t & f
    @test @inferred(-(f, t)) === StaticInt(f) - StaticInt(t)
    @test @inferred(+(f, t)) === StaticInt(f) + StaticInt(t)

    @test @inferred(^(t, f)) == ^(true, false)
    @test @inferred(^(t, t)) == ^(true, true)

    @test @inferred(^(2, f)) == 1
    @test @inferred(^(2, t)) == 2

    @test @inferred(^(BigInt(2), f)) == 1
    @test @inferred(^(BigInt(2), t)) == 2

    @test @inferred(div(t, t)) === t
    @test_throws DivideError div(t, f)

    @test @inferred(rem(t, t)) === f
    @test_throws DivideError rem(t, f)
    @test @inferred(mod(t, t)) === f

    @test @inferred(all((t, t, t)))
    @test !@inferred(all((t, f, t)))
    @test !@inferred(all((f, f, f)))

    @test @inferred(any((t, t, t)))
    @test @inferred(any((t, f, t)))
    @test !@inferred(any((f, f, f)))

    @test static(true) == true
    @test static(false) != true
    @test static(true) != false
    @test static(false) == false
    @test !(static(true) != true)
    @test !(static(false) == true)
    @test !(static(true) == false)
    @test !(static(false) != false)
    @test true == static(true)
    @test false != static(true)
    @test true != static(false)
    @test false == static(false)
    @test !(true != static(true))
    @test !(false == static(true))
    @test !(true == static(false))
    @test !(false != static(false))
end

@testset "operators" begin
    f = static(false)
    t = static(true)
    x = StaticInt(1)
    y = StaticInt(0)
    z = StaticInt(-1)
    @test @inferred(Static.eq(y)(x)) === f
    @test @inferred(Static.eq(x, x)) === t

    @test @inferred(Static.ne(y)(x)) === t
    @test @inferred(Static.ne(x, x)) === f

    @test @inferred(Static.gt(y)(x)) === t
    @test @inferred(Static.gt(y, x)) === f

    @test @inferred(Static.ge(y)(x)) === t
    @test @inferred(Static.ge(y, x)) === f

    @test @inferred(Static.lt(x)(y)) === t
    @test @inferred(Static.lt(x, y)) === f

    @test @inferred(Static.le(x)(y)) === t
    @test @inferred(Static.le(x, y)) === f

    @test @inferred(Static.ifelse(t, x, y)) === x
    @test @inferred(Static.ifelse(f, x, y)) === y

    sa1 = Static.add(x)
    sa0 = Static.add(y)
    saz = Static.add(z)
    da = Static.add(dynamic(x))

    sm1 = Static.mul(x)
    sm0 = Static.mul(y)
    smz = Static.mul(z)
    dm = Static.mul(dynamic(x))

    @test (sa1 ∘ sa1) === Static.add(static(2))
    @test (sm1 ∘ sa0) === sm1
    @test (sm0 ∘ sa0) === sm0
    @test (smz ∘ sa0) === smz
    @test (dm ∘ sa0) === dm
    @test (sm1 ∘ saz) === saz
    @test (sm0 ∘ saz) === sm0
    @test (sm1 ∘ da) === da
    @test (sm0 ∘ da) === sm0
    @test (sm1 ∘ smz) === smz  # 1 * -1
    @test (sm0 ∘ dm) === sm0
    @test (dm ∘ sm0) === sm0
    @test (sm1 ∘ dm) === dm
    @test (dm ∘ sm1) === dm
end

@testset "static interface" begin
    v = Val((:a, 1, true))

    @test static(1) === StaticInt(1)
    @test static(true) === True()
    @test static(:a) === StaticSymbol{:a}() === static("a") === static('a')
    @test Symbol(static(:a)) === :a
    @test static((:a, 1, true)) === (static(:a), static(1), static(true))
    @test @inferred(static(v)) === (static(:a), static(1), static(true))

    @test @inferred(Static.is_static(v)) === True()
    @test @inferred(Static.is_static(typeof(v))) === True()
    @test @inferred(Static.is_static(typeof(static(true)))) === True()
    @test @inferred(Static.is_static(typeof(static(1)))) === True()
    @test @inferred(Static.is_static(typeof(static(:x)))) === True()
    @test @inferred(Static.is_static(typeof(1))) === False()
    @test @inferred(Static.is_static(typeof((static(:x), static(:x))))) === True()
    @test @inferred(Static.is_static(typeof((static(:x), :x)))) === False()
    @test @inferred(Static.is_static(typeof(static(1.0)))) === True()

    @test @inferred(Static.known(v)) === (:a, 1, true)
    @test @inferred(Static.known(typeof(v))) === (:a, 1, true)
    @test @inferred(Static.known(typeof(static(true))))
    @test @inferred(Static.known(typeof(static(false)))) === false
    @test @inferred(Static.known(typeof(static(1.0)))) === 1.0
    @test @inferred(Static.known(typeof(static(1)))) === 1
    @test @inferred(Static.known(typeof(static(:x)))) === :x
    @test @inferred(Static.known(typeof(1))) === nothing
    @test @inferred(Static.known(typeof((static(:x), static(:x))))) === (:x, :x)
    @test @inferred(Static.known(typeof((static(:x), :x)))) === (:x, nothing)

    @test @inferred(Static.dynamic((static(:a), static(1), true))) === (:a, 1, true)
end

@testset "promote_shape" begin
    x = (static(1), 1)
    y = (1, static(1), 1)
    @test @inferred(Base.promote_shape(x, x)) === (static(1), 1)
    @test @inferred(Base.promote_shape(x, y)) === (static(1), static(1), static(1))
    @test @inferred(Base.promote_shape(y, x)) === (static(1), static(1), static(1))
    @test static_promote(1, nothing) === 1
    @test static_promote(nothing, 1) === 1
    @test static_promote(nothing, nothing) === nothing
    @test_throws ErrorException Base.promote_shape((static(1),), (static(2),))
end

@testset "tuple utilities" begin
    x = (static(1), static(2), static(3))
    y = (static(3), static(2), static(1))
    z = (static(1), static(2), static(3), static(4))
    T = Tuple{Int, Float64, String}
    @test @inferred(Static.invariant_permutation(x, x)) === True()
    @test @inferred(Static.invariant_permutation(x, y)) === False()
    @test @inferred(Static.invariant_permutation(x, z)) === False()

    @test @inferred(Static.permute(x, Val(x))) === x
    @test @inferred(Static.permute(x, (static(1), static(2)))) === (static(1), static(2))
    @test @inferred(Static.permute(x, x)) === x
    @test @inferred(Static.permute(x, y)) === y
    @test @inferred(Static.eachop(getindex, x)) === x

    @test Static.field_type(typeof((x = 1, y = 2)), :x) <: Int
    @test Static.field_type(typeof((x = 1, y = 2)), static(:x)) <: Int
    function get_tuple_add(::Type{T}, ::Type{X}, dim::StaticInt) where {T, X}
        Tuple{Static.field_type(T, dim), X}
    end
    @test @inferred(Static.eachop_tuple(Static.field_type, y, T)) ===
          Tuple{String, Float64, Int}
    @test @inferred(Static.eachop_tuple(get_tuple_add, y, T, String)) ===
          Tuple{Tuple{String, String}, Tuple{Float64, String}, Tuple{Int, String}}
    @test @inferred(Static.find_first_eq(static(1), y)) === static(3)
    @test @inferred(Static.find_first_eq(static(2), (1, static(2)))) === static(2)
    @test @inferred(Static.find_first_eq(static(2), (1, static(2), 3))) === static(2)
    # inferred is Union{Int,Nothing}
    @test Static.find_first_eq(1, map(Int, y)) === 3

    @testset "reduce_tup" begin
        for n in 2:16
            x = ntuple(
                _ -> rand(Bool) ? rand() : (rand(Bool) ? rand(0x00:0x1f) : rand(0:31)),
                n)
            @test @inferred(Static.reduce_tup(+, x)) ≈ reduce(+, x)
        end
    end
end

@testset "invperm" begin
    perm = static((10, 3, 4, 5, 6, 2, 9, 7, 8, 1))
    invp = static((10, 6, 2, 3, 4, 5, 8, 9, 7, 1))
    @test @inferred(invperm(perm)) === invp
end

@testset "NDIndex" begin
    x = NDIndex((1, 2, 3))
    y = NDIndex((1, static(2), 3))
    z = NDIndex(static(3), static(3), static(3))

    @testset "constructors" begin
        @test static(CartesianIndex(3, 3, 3)) === z ==
              Base.setindex(Base.setindex(x, 3, 1), 3, 2)
        @test @inferred(CartesianIndex(z)) === @inferred(Static.dynamic(z)) ===
              CartesianIndex(3, 3, 3)
        @test @inferred(Static.known(z)) === (3, 3, 3)
        @test Tuple(@inferred(NDIndex{0}())) === ()
        @test @inferred(NDIndex{3}(1, static(2), 3)) === y
        @test @inferred(NDIndex{3}((1, static(2), 3))) === y
        @test @inferred(NDIndex{3}((1, static(2)))) === NDIndex(1, static(2), static(1))
        @test @inferred(NDIndex(x, y)) === NDIndex(1, 2, 3, 1, static(2), 3)
    end

    @test @inferred(Base.IteratorsMD.split(x, Val(2))) === (NDIndex(1, 2), NDIndex(3))
    @test @inferred(length(x)) === 3
    @test @inferred(length(typeof(x))) === 3
    NDIndex2{I <: Tuple{Vararg{Union{StaticInt, Int}, 2}}} = NDIndex{2, I}
    @test length(NDIndex2) === 2
    @test @inferred(y[2]) === 2
    @test @inferred(y[static(2)]) === static(2)

    @test @inferred(-y) === NDIndex((-1, -static(2), -3))
    @test @inferred(y+y) === NDIndex((2, static(4), 6))
    @test @inferred(y-y) === NDIndex((0, static(0), 0))
    @test @inferred(zero(x)) === NDIndex(static(0), static(0), static(0))
    @test @inferred(oneunit(x)) === NDIndex(static(1), static(1), static(1))
    @test @inferred(x*3) === NDIndex((3, 6, 9))
    @test @inferred(3*x) === NDIndex((3, 6, 9))

    @test @inferred(min(x, z)) === x
    @test @inferred(max(x, z)) === NDIndex(3, 3, 3)
    @test !@inferred(isless(y, x))
    @test @inferred(isless(x, z))
    @test @inferred(Static.lt(oneunit(z), z)) === static(true)

    A = rand(3, 3, 3)
    @test @inferred(to_indices(A, axes(A), (x,))) === (1, 2, 3)
    @test @inferred(to_indices(A, axes(A), ([y, y],))) == ([y, y],)

    # issue #44
    @test deleteat!(Union{}[], Union{}[]) == Union{}[]
end

# for some reason this can't be inferred when in the "Static.jl" test set
known_length(x) = known_length(typeof(x))
known_length(::Type{T}) where {N, T <: Tuple{Vararg{Any, N}}} = N
known_length(::Type{T}) where {T} = nothing
maybe_static_length(x) = Static.maybe_static(known_length, length, x)
x = ntuple(+, 10)
y = 1:10
@test @inferred(maybe_static_length(x)) === StaticInt(10)
@test @inferred(maybe_static_length(y)) === 10

@testset "StaticFloat64" begin
    f = static(1.0)
    @test @inferred(dynamic(f)) === @inferred(known(f)) === 1.0
    for i in -10:10
        for j in -10:10
            @test i + j == @inferred(Static.StaticInt(i)+Static.StaticFloat64(j)) ==
                  @inferred(i+Static.StaticFloat64(Static.StaticFloat64(j))) ==
                  @inferred(Static.StaticFloat64(i)+j) ==
                  @inferred(Static.StaticFloat64(i)+Static.StaticInt(j)) ==
                  @inferred(Static.StaticFloat64(i)+Static.StaticFloat64(j))
            @test i - j == @inferred(Static.StaticInt(i)-Static.StaticFloat64(j)) ==
                  @inferred(i-Static.StaticFloat64(j)) ==
                  @inferred(Static.StaticFloat64(i)-Static.StaticInt(j)) ==
                  @inferred(Static.StaticFloat64(i)-j) ==
                  @inferred(Static.StaticFloat64(i)-Static.StaticFloat64(j))
            @test i * j == @inferred(Static.StaticInt(i)*Static.StaticFloat64(j)) ==
                  @inferred(i*Static.StaticFloat64(j)) ==
                  @inferred(Static.StaticFloat64(i)*Static.StaticInt(j)) ==
                  @inferred(Static.StaticFloat64(i)*j) ==
                  @inferred(Static.StaticFloat64(i)*Static.StaticFloat64(j))
            i == j == 0 && continue
            @test i / j == @inferred(Static.StaticInt(i)/Static.StaticFloat64(j)) ==
                  @inferred(i/Static.StaticFloat64(j)) ==
                  @inferred(Static.StaticFloat64(i)/Static.StaticInt(j)) ==
                  @inferred(Static.StaticFloat64(i)/j) ==
                  @inferred(Static.StaticFloat64(i)/Static.StaticFloat64(j))
        end
        if i ≥ 0
            @test sqrt(i) == @inferred(sqrt(Static.StaticInt(i))) ==
                  @inferred(sqrt(Static.StaticFloat64(i))) ==
                  @inferred(sqrt(Static.StaticFloat64(Float64(i))))
        end
    end
    @test Static.floortostaticint(1.0) === 1
    @test Static.floortostaticint(prevfloat(2.0)) === 1
    @test @inferred(Static.floortostaticint(Static.StaticFloat64(1.0))) ===
          Static.StaticInt(1)
    @test @inferred(Static.floortostaticint(Static.StaticFloat64(prevfloat(2.0)))) ===
          Static.StaticInt(1)

    @test Static.roundtostaticint(1.0) === 1
    @test Static.roundtostaticint(prevfloat(2.0)) === 2
    @test @inferred(Static.roundtostaticint(Static.StaticFloat64(1.0))) ===
          Static.StaticInt(1)
    @test @inferred(Static.roundtostaticint(Static.StaticFloat64(prevfloat(2.0)))) ===
          Static.StaticInt(2)
    @test @inferred(round(Static.StaticFloat64{1.0}())) === Static.StaticFloat64(1)
    @test @inferred(round(Static.StaticFloat64(prevfloat(2.0)))) ===
          Static.StaticFloat64(ComplexF64(2))

    x = static(5.0)
    y = static(10.0)
    @test @inferred(rem(x, y)) === x
    @test @inferred(rem(x, 10.0)) === 5.0
    @test @inferred(rem(5.0, y)) === 5.0

    @test @inferred(min(x, y)) === x
    @test @inferred(min(x, 10.0)) === 5.0
    @test @inferred(min(5.0, y)) === 5.0

    @test @inferred(max(x, y)) === y
    @test @inferred(max(x, 10.0)) === 10.0
    @test @inferred(max(5.0, y)) === 10.0

    @test @inferred(isless(x, y))
    @test @inferred(isless(x, 10.0))
    @test @inferred(isless(5.0, y))

    fone = static(1.0)
    fzero = static(0.0)
    @test @inferred(isone(fone))
    @test @inferred(isone(one(fzero)))
    @test @inferred(isone(fzero)) === false

    @test @inferred(iszero(fone)) === false
    @test @inferred(iszero(fzero))
    @test @inferred(iszero(zero(typeof(fzero))))

    @test typeof(fone)(1) isa Static.StaticFloat64
    @test typeof(fone)(1.0) isa Static.StaticFloat64

    @test @inferred(eltype(Static.StaticFloat64(static(1)))) <: Float64
    @test @inferred(promote_rule(typeof(fone), Int)) <: promote_type(Float64, Int)
    @test @inferred(promote_rule(typeof(fone), Float64)) <: Float64
    @test @inferred(promote_rule(typeof(fone), Float32)) <: Float32
    @test @inferred(promote_rule(typeof(fone), Float16)) <: Float16

    @test @inferred(inv(static(2.0))) === static(inv(2.0)) === inv(static(2))

    @test @inferred(static(2.0)^2.0) === 2.0^2.0

    @testset "trig" begin
        for f in [sin, cos, tan, asin, atan, acos, sinh, cosh, tanh,
            asinh, atanh, exp, exp2,
            exp10, expm1, log, log2, log10, log1p, exponent, sqrt, cbrt, sec, csc, cot, sech,
            secd, csch, cscd, cotd, cosd, tand, asind, acosd, atand, acotd, sech, coth, asech,
            acsch, deg2rad, mod2pi, sinpi, cospi]
            @info "Testing $f(0.5)"
            @test @inferred(f(static(0.5))) === static(f(0.5))
        end
    end
    @test @inferred(asec(static(2.0))) === static(asec(2.0))
    @test @inferred(acsc(static(2.0))) === static(acsc(2.0))
    @test @inferred(acot(static(2.0))) === static(acot(2.0))
    @test @inferred(asecd(static(2.0))) === static(asecd(2.0))
    @test @inferred(acscd(static(2.0))) === static(acscd(2.0))
    @test @inferred(acoth(static(2.0))) === static(acoth(2.0))
    @info "Testing acosh(1.5)"
    @inferred acosh(static(1.5))

    @test @inferred(isinteger(static(1.0)))
    @test @inferred(!isinteger(static(1.5)))
end

@testset "string/print/show" begin
    f = static(float(2))
    repr(f)
    @test repr(static(float(1))) == "static($(float(1)))"
    @test repr(static(1)) == "static(1)"
    @test repr(static(:x)) == "static(:x)"
    @test repr(static(true)) == "static(true)"
    @test repr(static(CartesianIndex(1, 1))) == "NDIndex(static(1), static(1))"
    @test string(static(true)) == "static(true)" == "$(static(true))"
    @test repr(static(1):static(10)) == "static(1):static(10)"
    @test repr(static(1):static(2):static(9)) == "static(1):static(2):static(9)"
end

include("ranges.jl")

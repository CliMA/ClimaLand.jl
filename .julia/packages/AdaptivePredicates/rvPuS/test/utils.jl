abstract type Method end
struct MacroMethod <: Method
    f::Symbol
    nargs::Int
end
struct ArithmeticMethod <: Method
    f::Symbol
    nargs::Int
    argdims::Tuple
    # argdims:
    #   d > 0 ⟹ NTuple{d, Float}. If d == 1, just becomes a Float
    #   d < 0 ⟹ random vector of length rand(1:8)
    #   d isa Function ⟹ A function of args such that d(args) returns the dimension, and the result returned if zeros(T, d(args))
    ArithmeticMethod(f, nargs, argdims::Tuple) = (length(argdims) == nargs || throw(ArgumentError("length(argdims) ≠ nargs"))) && return new(f, nargs, argdims)
end
struct PredicateMethod <: Method
    f::Symbol
    nargs::Int
    argdims::Int
end
argdims(method::MacroMethod, i) = 1
argdims(method::ArithmeticMethod, i) = method.argdims[i]
argdims(method::PredicateMethod, i) = method.argdims

function _rand(::Type{T}) where {T}
    if rand() < 1 / 4
        # Generate a float
        s = rand((-1, 1))
        f = rand(T)
        e = T == Float64 ? rand(-100:200) : rand(-24:24)
        return s * (one(f) + f) * T(2)^e # no subnormal normals will be generated here
    else
        # Generate an integer 
        return T(rand(-100:100))
    end
end
const NTESTS = 50000
macro repeat(args...)
    if length(args) == 1
        n = NTESTS
        expr = args[1]
    elseif length(args) == 2
        n, expr = args[1:2]
    else
        throw(ArgumentError("Invalid expression passed to @repeat."))
    end
    quote
        foreach(_ -> $(esc(expr)), 1:$(esc(n)))
    end
end

function genargs(method::Method, ::Type{T}) where {T}
    args = ()
    for i in 1:method.nargs
        argdim = argdims(method, i)
        arg = if argdim == 1
            _rand(T)
        elseif argdim isa Function
            zeros(T, argdim(args))
        elseif argdim > 0
            ntuple(_ -> _rand(T), Val(argdim))
        else
            [_rand(T) for _ in 1:rand(1:8)]
        end
        args = (args..., arg)
    end
    if method isa PredicateMethod
        f = String(method.f)
        if endswith(f, "adapt")
            args = (args..., _rand(T))
        end
        if f ∈ ("incircle", "orient2") && rand() < 1 / 2 # 50% chance we test Complex instead of NTuple{2}
            new_args = ()
            for arg in args
                if arg isa NTuple{2}
                    new_args = (new_args..., Complex(arg[1], arg[2]))
                else
                    new_args = (new_args..., arg)
                end
            end
            args = new_args
        end
    end
    return add_lengths(method, args)
end
function add_lengths(method, args) # adds elen, flen 
    new_args = ()
    for i in eachindex(args)
        if args[i] isa Vector && (i ≠ lastindex(args) || method.f == :estimate)
            new_args = (new_args..., length(args[i]), args[i])
        else
            new_args = (new_args..., args[i])
        end
    end
    return new_args
end
function _retup_complex(args)
    new_args = ()
    for arg in args
        if arg isa Complex
            new_args = (new_args..., (arg.re, arg.im))
        else
            new_args = (new_args..., arg)
        end
    end
    return new_args
end

struct Expansion{T}
    h::Vector{T}
    hlen::Int # not necessarily length(h)
    function Expansion(h::Vector{T}, hlen::Int; check_zero=false) where {T}
        check_zero && @test all(iszero, view(h, hlen+1:length(h))) && hlen ≤ length(h) # For compress, the remaining components don't need to be zero
        return new{T}(h, hlen)
    end
end
Expansion(x::Number; kwargs...) = Expansion([x], 1; kwargs...)
Expansion(x::NTuple{N,<:AbstractFloat}; kwargs...) where {N} = Expansion(collect(x), N; kwargs...)
Expansion(xh::Tuple{NTuple{N,<:AbstractFloat},Int}; kwargs...) where {N} = Expansion(collect(xh[1]), xh[2]; kwargs...)
Expansion(xh::Tuple{<:Vector,Int}; kwargs...) = Expansion(xh[1], xh[2]; kwargs...)

expand(h::Expansion) = view(h.h, 1:h.hlen)
Base.sum(h::Expansion) = sum(expand(h))
function Base.:(==)(e::Expansion, f::Expansion)
    return sum(e) == sum(f)
end

function ⪧(LHS, RHS; kwargs...)
    return Expansion(LHS; kwargs...) == Expansion(RHS; kwargs...)
end

struct FailedTest
    LHSf::Function
    RHSf::Function
    LHSv::Any
    RHSv::Any
    args::Tuple
end

function retry(failures::Vector{FailedTest})
    remaining = Vector{FailedTest}()
    foreach(failures) do failure 
        LHSf = failure.LHSf
        RHSv, args = failure.RHSv, failure.args
        LHSv = LHSf(args...)
        flag = ⪧(LHSv, RHSv; check_zero = nameof(LHSf) ∉ non_check_zeros)
        if !flag
            push!(remaining, failure)
        end
    end
    return remaining
end

const non_check_zeros = (:compress, :expansion_sum_zeroelim2, :expansion_sum_zeroelim1)
function test_f(method::Method; failures)
    return _test_f(method; failures)
end
function test_f(method::PredicateMethod; failures)
    for suffix in (:fast, :exact, :adapt, :slow, Symbol())
        g = Symbol(method.f, suffix)
        _method = PredicateMethod(g, method.nargs, method.argdims)
        _test_f(_method; failures)
    end
    return true
end
function _test_f(method::Method; failures)
    args64 = genargs(method, Float64)
    args32 = genargs(method, Float32)
    LHS = getproperty(AP, method.f)
    cmethod = method isa PredicateMethod ? Symbol(replace(String(method.f), "2" => "2d", "3" => "3d")) : method.f
    RHS = getproperty(C, cmethod)
    LHSf64 = @inferred LHS(args64...)
    LHSf32 = @inferred LHS(args32...)
    RHSf64 = RHS(_retup_complex(args64)...)
    RHSf32 = RHS(_retup_complex(args32)...)
    flag1 = ⪧(LHSf64, RHSf64; check_zero=method.f ∉ non_check_zeros)
    flag2 = ⪧(LHSf32, RHSf32; check_zero=method.f ∉ non_check_zeros)
    !flag1 && push!(failures, FailedTest(LHS, RHS, LHSf64, RHSf64, args64))
    !flag2 && push!(failures, FailedTest(LHS, RHS, LHSf32, RHSf32, args32))
    @test LHSf64 ⪧ RHSf64 check_zero = method.f ∉ non_check_zeros
    @test LHSf32 ⪧ RHSf32 check_zero = method.f ∉ non_check_zeros
    if method isa PredicateMethod
        if method.f ∈ (:orient2, :incircle, :orient3, :insphere)
            LHSp = getproperty(AP, Symbol(method.f, :p))
            LHSf64p::Int = @inferred LHSp(args64...)
            LHSf32p::Int = @inferred LHSp(args32...)
            flag1 = LHSf64p == Int(sign(RHSf64))
            flag2 = LHSf32p == Int(sign(RHSf32))
            !flag1 && push!(failures, FailedTest(LHSp, RHS, LHSf64p, Int(sign(RHSf64)), args64))
            !flag2 && push!(failures, FailedTest(LHSp, RHS, LHSf32p, Int(sign(RHSf32)), args32))
            @test LHSf64p == Int(sign(RHSf64))
            @test LHSf32p == Int(sign(RHSf32))
        end
        @test LHSf64 isa Float64
        @test LHSf32 isa Float32
        # Test concurrency
        if rand() < 100 / NTESTS
            ap64 = Vector{Float64}(undef, 24)
            ap32 = Vector{Float32}(undef, 24)
            Base.Threads.@threads for i in 1:24
                ap64[i] = LHS(args64...)
                ap32[i] = LHS(args32...)
            end
            @test all(==(RHSf64), ap64)
            @test all(==(RHSf32), ap32)
        end
    end
    return true
end
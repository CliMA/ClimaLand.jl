# SPDX-License-Identifier: MIT
# Author: Pierre Lairez

module Codegen

using IntervalArithmetic
using StaticArrays

import Base: +, *, -, one, convert, promote_rule

export coord, @genpredicate

@enum FilterEnum fastfp_flt zerotest_flt accuratefp_flt interval_flt exact_flt naive_flt nothing_flt

struct InvalidPrecisionError{T} <: Exception
    x::T
end
function Base.showerror(io :: IO, e :: InvalidPrecisionError{T}) where {T}
    print(io, "Invalid precision used for input ", e.x, "::$T. You can only use Float64 numbers in predicates.")
end

@inline _coord(x :: NTuple{N,Float64}) where {N} = x
_coord(x :: Tuple) = throw(InvalidPrecisionError(x))
_coord(x) = coord(x) 

@inline coord(x) = _coord(Tuple(x))
@inline coord(x :: Complex) = _coord(reim(x))

function adddict(a, b)
    d = copy(a)
    for k in keys(b)
        if haskey(a, k)
            d[k] += b[k]
        else
            d[k] = b[k]
        end
    end
    return d
end

function mergedict(a, b)
    d = copy(a)
    for k in keys(b)
        if haskey(a, k)
            d[k] = union(d[k], b[k])
        else
            d[k] = b[k]
        end
    end
    return d
end

mutable struct Formula <: Number
    id :: Symbol
    head :: Symbol
    args :: Vector
    group :: Union{Nothing, Symbol}
end

function Formula(name :: Union{Expr, Symbol}, group = nothing)
    s = gensym()
    Formula(s, :sym, [name], group)
end

macro var(args...)
    vars = [Formula(v) for v in args]
    Expr(:block,
         [Expr(:(=), esc(args[i]), vars[i]) for i in 1:length(vars)]...)
end

macro point2(args...)
    vars = [SVector(Formula(:($v[1])), Formula(:($v[2]))) for v in args]
    Expr(:block,
         [Expr(:(=), esc(args[i]), vars[i]) for i in 1:length(vars)]...)
end

function +(f :: Formula, g :: Formula)
    Formula(gensym(), :+, [f, g], nothing)
end

function *(f :: Formula, g :: Formula)
    Formula(gensym(), :*, [f, g], nothing)
end

function -(f :: Formula, g :: Formula)
    Formula(gensym(), :-, [f, g], nothing)
end

function -(f :: Formula)
    Formula(gensym(), :-, [f], nothing)
end



function convert(::Type{Formula}, i :: T) where T <: Integer
    Formula(gensym(), :const, [i], nothing)
end

zero(::Type{Formula}) = convert(Formula, 0)

promote_rule(::Type{T}, ::Type{Formula}) where {T <: Integer} = Formula
promote_rule(::Type{Bool}, ::Type{Formula}) = Formula # method ambiguity


function group!(args...)
    @gensym g
    for f in args
        f.group = g
    end
end


function accumulator(f :: Formula)
    if f.head == :sym
        ret = (deg = Dict(), groups = Dict(), bound = Inf, error = 0.0)
    elseif f.head == :+ || (f.head == :- && length(f.args) == 2)
        af = accumulator(f.args[1])
        ag = accumulator(f.args[2])
        @assert af.deg == ag.deg

        bound = nextfloat(af.bound + ag.bound)
        error = nextfloat(nextfloat(af.error + ag.error) + eps(bound)/2)
        ret = (deg = af.deg, groups = mergedict(af.groups, ag.groups), bound = bound, error = error)
    elseif f.head == :*
        af = accumulator(f.args[1])
        ag = accumulator(f.args[2])
        bound = nextfloat(af.bound * ag.bound)
        error = nextfloat(nextfloat(nextfloat(af.error*ag.bound) + nextfloat(af.bound*ag.error)) + eps(bound)/2)
        ret = (deg = adddict(af.deg, ag.deg), groups = mergedict(af.groups, ag.groups), bound = bound, error = error)
    elseif f.head == :-
        ret = accumulator(f.args[1])
    elseif f.head == :const
        ret = (deg = Dict(), groups = Dict(), bound = abs(float(f.args[1])), error = 0.0)
    else
        throw(DomainError())
    end

    if !isnothing(f.group)
        @assert f.head == :sym || f.head == :- && f.args[1].head == f.args[2].head == :sym
        ret = (deg = Dict(f.group => 1), groups = Dict(f.group => Set([f.id])), bound = 1.0, error = eps(1.0)/2)
    end

    return ret
end


function evalcode(f :: Formula, conv = nothing)
    code = []
    stack = [f]
    syms = Set{Symbol}()

    while !isempty(stack)
        e = pop!(stack)

        if e.id ∈ syms
            continue
        end

        s = (e.id)

        if e.head == :sym
            if isnothing(conv)
                push!(code, Expr(:(=), s, e.args[1]))
            else
                push!(code, Expr(:(=), s, conv(e.args[1])))
            end
        elseif e.head == :const
            push!(code, code, Expr(:(=), s, e.args[1]))
        else
            if all(a.id ∈ syms for a in e.args)
                push!(code, Expr(:(=), s, Expr(:call, e.head, ((a.id) for a in e.args)...)))
            else
                push!(stack, e, e.args...)
                continue
            end
        end

        push!(syms, e.id)
    end

    return quote
        $(code...)
        $(f.id)
    end
end


function fastfilter(f :: Formula ; withretcode :: Bool = false)
    code = []

    @gensym res

    let
        fpcode = quote
            $res = $(evalcode(f))
        end
        push!(code, fpcode)
    end

    # Now code contains the computation of the formula, nothing more. We can
    # reliably determine the sign of the result if its absolute value is larger
    # than acc.error*scaling. The scaling depends on the homogeneity structure
    # of the formula.

    # The format of acc.groups is Dict(g1=>[a,b,c,...], g2=>[d,e,f,...], ...),
    # where a, b, c, ... are input gates (or difference or input gates, which
    # are considered as new input gates). The formula is expected to be
    # multihomogeneous w.r.t. to each group. The degree of homogeneity w.r.t to
    # the group g is acc.deg[g]. Therefore, the appropriate scaling is the
    # product of all λ[g]^acc.deg[g], where λ[g] is the maximum
    # absolute value of an input in group g.

    # The accurate way of determining if abs(res) >
    # acc.error*Π λ[g]^acc.deg[g] is to do the computation with fp numbers.
    # But we can check quickly a sufficient condition just with integer
    # arithmetic on the exponents.

    # There are two fundamental things to check: that no overflow occurred with computing res
    # and that no underflow occurs when computing acc.error*Π λ[g]^acc.deg[g].

    acc = accumulator(f)
    totdeg = UInt(sum(values(acc.deg)))

    @gensym signres

    let
        # we first try to determine if the absolute value of the result is
        # bigger than the maximal possible absolute error by looking only at the
        # exponents.

        emask = Base.exponent_mask(Float64)
        mantissabits = UInt(Base.Math.significand_bits(Float64))
        grouplogs = []
        for (idx, g) in acc.groups
            @gensym glog
            push!(grouplogs, glog)

            alllogs = [:(reinterpret(UInt64, $v) & $emask) for v in g]
            push!(code, :($glog = max($(alllogs...)) >> $mantissabits))
            # we could probably work without the shift, but it would make it
            # harder to think about overflows/underflows.

            if acc.deg[idx] != 1
                push!(code, :($glog *= $(acc.deg[idx])))
            end
        end

        errlog = ((reinterpret(UInt64, acc.error) & emask) >> mantissabits) -
            totdeg * ( Base.exponent_one(Float64) >> mantissabits ) + totdeg

        @gensym ε
        @gensym reslog
        @gensym rawres

        filter = quote
            $ε = $(Expr(:call, :+, errlog, grouplogs...))
            $rawres = reinterpret(UInt64, $res)
            $signres = Base.ifelse($rawres & $(Base.sign_mask(Float64)) == 0, 1, -1)
            $reslog = ($rawres & $emask) >> $mantissabits

            if $reslog != $(emask >> mantissabits) && # protect against Inf and Nan
                #$ε ≥ 0 &&       # protect against underflow in computation of ε
                $reslog > $ε    # main test
                # we should check that ε is not
                # negative. But since we are in unsigned arithmetic, ε < 0 means
                # actually that ε is big, so reslog > ε cannot hold.

                $(withretcode ? :(return ($signres, $fastfp_flt)) : :(return $signres))
            end
        end
        push!(code, filter)
    end


    let
        # Quick test is non conclusive. We check if there is a group in which
        # all the variables are zero, in which case the result is zero.

        filter = quote
            if $(Expr(:call, :|,
                      (Expr(:call, :&,
                            (Expr(:call, :(==), v, 0.0)
                             for v in g)...)
                       for g in values(acc.groups))...)
                 )

                $(withretcode ? :(return (0, $zerotest_flt)) : :(return 0))
            end
        end
        push!(code, filter)
    end


    let
        # We now go to the finer, more straight forward, test.

        groupabs = []
        for (idx, g) in acc.groups
            @gensym gabs
            push!(groupabs, gabs)
            allabs = [:(abs($v)) for v in g]
            push!(code, :($gabs = max($(allabs...))))

            if acc.deg[idx] != 1
                push!(code, :($gabs = $gabs^$(acc.deg[idx])))
            end
        end

        errabs =  acc.error * (1+eps(1.0))^totdeg

        @gensym ε

        filter = quote
            $ε = *($(groupabs...)) * $errabs
            if !issubnormal($ε) && $ε > 0 && isfinite($res) && abs($res) > $ε
                $(withretcode ? :(return ($signres, $accuratefp_flt)) : :(return $signres))
            end
        end
        push!(code, filter)
    end

    return quote
        $(code...)
    end

end


function ivfilter(f :: Formula ; withretcode :: Bool = false)

    @gensym ivres
    quote
        # We now resort to interval arithmetic It is an interesting filter when
        # the data is made of exactly representable integers.
        $ivres = $(evalcode(f,  s -> :( interval($s) )))
        if isstrictless($ivres, interval(0.0))
            $(withretcode ? :(return (-1, $interval_flt)) : :(return -1))
        elseif isstrictless(interval(0.0), $ivres)
            $(withretcode ? :(return (1, $interval_flt)) : :(return 1))
        elseif mag($ivres) == 0
            $(withretcode ? :(return (0, $interval_flt)) : :(return 0))
        end
    end

end

function exfilter(f :: Formula ; withretcode :: Bool = false)

    @gensym exrec
    quote
        # Exact arithmetic. Always conclusive.
        $exrec = $(evalcode(f, s -> :( Rational{BigInt}($s) )))
        $(withretcode ? :(return (Int(sign($exrec)), $exact_flt)) : :(return Int(sign($exrec))))
    end

end

function naivefilter(f :: Formula  ; withretcode :: Bool = false)
    @gensym fpres
    quote
        # Exact arithmetic. Always conclusive.
        $fpres = $(evalcode(f))
        if $fpres < 0
            $(withretcode ? :(return (-1, $naive_flt)) : :(return -1))
        elseif $fpres > 0
            $(withretcode ? :(return (1, $naive_flt)) : :(return 1))
        else
            $(withretcode ? :(return (0, $naive_flt)) : :(return 0))
        end
    end
end

"""

Generate sign predicate for a function that computes a polynomial in the
coordinates of the arguments.

"""
macro genpredicate(args...)
    if first(args) == :nogeneric
        defgeneric = false
        fun = args[2]
    else
        defgeneric = true
        fun = args[1]
    end

    sig = fun.args[1]
    args = sig.args[2:length(sig.args)]

    @assert fun.head == :function

    nargs = []
    tupleconv = []
    input = []
    for a in args
        if isa(a, Symbol)
            push!(input, Formula(a))
            push!(nargs, :($a :: Float64))
        elseif isa(a, Expr) && a.head == :(::)
            v, dim = a.args
            push!(input, SVector((Formula(:($v[$i])) for i in 1:dim)...))
            push!(nargs, :($v :: NTuple{$dim, Float64}))
            push!(tupleconv, :($v = _coord($v)))
        else
            throw(DomainError("Unknown argument $a"))
        end
    end

    base = string(sig.args[1])
    mainf = esc(Symbol(base))
    slowf = esc(Symbol(base, "_slow"))
    referencef = esc(Symbol(base, "_reference"))
    naivef = esc(Symbol(base, "_naive"))

    debug = s -> esc(Symbol(s.args[1], "_dbg"))

    nsig = (a for a in nargs)
    vars = (a.args[1] for a in nargs)

    # This is the fun part: we inject into the body of the function given in
    # argument the formal input. This gives an object of type Formula that
    # represents the polynomial.
    formula = Core.eval(__module__,
         Expr(:call,
              Expr(:->, Expr(:tuple, vars...), fun.args[2]),
              input...))

    if defgeneric
        genfun = quote
            function $(mainf)($(vars...))
                $(tupleconv...)
                return $(mainf)($(vars...))
            end
        end
    else
        genfun = quote end
    end


    quote
        function $(naivef)($(nsig...))
            $(naivefilter(formula))
        end

        function $(debug(naivef))($(nsig...))
            $(naivefilter(formula, withretcode=true))
        end

        function $(referencef)($(nsig...))
            $(exfilter(formula))
        end

        function $(debug(referencef))($(nsig...))
            $(exfilter(formula, withretcode=true))
        end

        function $(slowf)($(nsig...))
            $(ivfilter(formula))
            return $(referencef)($(vars...))
        end

        function $(debug(slowf))($(nsig...))
            $(ivfilter(formula, withretcode=true))
            return $(debug(referencef))($(vars...))
        end

        Core.@__doc__(@inline function $(mainf)($(nsig...))
            $(fastfilter(formula))
            return $(slowf)($(vars...))
        end)

        function $(debug(mainf))($(nsig...))
            $(fastfilter(formula, withretcode=true))
            return $(debug(slowf))($(vars...))
        end

        $(genfun)

    end




end




end

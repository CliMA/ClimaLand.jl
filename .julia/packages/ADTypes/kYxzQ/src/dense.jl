"""
    AutoChainRules{RC}

Struct used to select an automatic differentiation backend based on [ChainRulesCore.jl](https://github.com/JuliaDiff/ChainRulesCore.jl) (see the list [here](https://juliadiff.org/ChainRulesCore.jl/stable/index.html#ChainRules-roll-out-status)).

Defined by [ADTypes.jl](https://github.com/SciML/ADTypes.jl).

# Constructors

    AutoChainRules(; ruleconfig)

# Fields

  - `ruleconfig::RC`: a [`ChainRulesCore.RuleConfig`](https://juliadiff.org/ChainRulesCore.jl/stable/rule_author/superpowers/ruleconfig.html) object.
"""
Base.@kwdef struct AutoChainRules{RC} <: AbstractADType
    ruleconfig::RC
end

mode(::AutoChainRules) = ForwardOrReverseMode()  # specialized in the extension

function Base.show(io::IO, backend::AutoChainRules)
    print(io, AutoChainRules, "(ruleconfig=", repr(backend.ruleconfig; context = io), ")")
end

"""
    AutoDiffractor

Struct used to select the [Diffractor.jl](https://github.com/JuliaDiff/Diffractor.jl) backend for automatic differentiation.

Defined by [ADTypes.jl](https://github.com/SciML/ADTypes.jl).

# Constructors

    AutoDiffractor()
"""
struct AutoDiffractor <: AbstractADType end

mode(::AutoDiffractor) = ForwardOrReverseMode()

"""
    AutoEnzyme{M,A}

Struct used to select the [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl) backend for automatic differentiation.

Defined by [ADTypes.jl](https://github.com/SciML/ADTypes.jl).

# Constructors

    AutoEnzyme(; mode::M=nothing, function_annotation::Type{A}=Nothing)

# Type parameters

  - `A` determines how the function `f` to differentiate is passed to Enzyme. It can be:

      + a subtype of `EnzymeCore.Annotation` (like `EnzymeCore.Const` or `EnzymeCore.Duplicated`) to enforce a given annotation
      + `Nothing` to simply pass `f` and let Enzyme choose the most appropriate annotation

# Fields

  - `mode::M` determines the autodiff mode (forward or reverse). It can be:

      + an object subtyping `EnzymeCore.Mode` (like `EnzymeCore.Forward` or `EnzymeCore.Reverse`) if a specific mode is required
      + `nothing` to choose the best mode automatically
"""
struct AutoEnzyme{M, A} <: AbstractADType
    mode::M
end

function AutoEnzyme(;
        mode::M = nothing, function_annotation::Type{A} = Nothing) where {M, A}
    return AutoEnzyme{M, A}(mode)
end

mode(::AutoEnzyme) = ForwardOrReverseMode()  # specialized in the extension

function Base.show(io::IO, backend::AutoEnzyme{M, A}) where {M, A}
    print(io, AutoEnzyme, "(")
    !isnothing(backend.mode) && print(io, "mode=", repr(backend.mode; context = io))
    !isnothing(backend.mode) && !(A <: Nothing) && print(io, ", ")
    !(A <: Nothing) && print(io, "function_annotation=", repr(A; context = io))
    print(io, ")")
end

"""
    AutoFastDifferentiation

Struct used to select the [FastDifferentiation.jl](https://github.com/brianguenter/FastDifferentiation.jl) backend for automatic differentiation.

Defined by [ADTypes.jl](https://github.com/SciML/ADTypes.jl).

# Constructors

    AutoFastDifferentiation()
"""
struct AutoFastDifferentiation <: AbstractADType end

mode(::AutoFastDifferentiation) = SymbolicMode()

"""
    AutoFiniteDiff{T1,T2,T3}

Struct used to select the [FiniteDiff.jl](https://github.com/JuliaDiff/FiniteDiff.jl) backend for automatic differentiation.

Defined by [ADTypes.jl](https://github.com/SciML/ADTypes.jl).

# Constructors

    AutoFiniteDiff(;
        fdtype=Val(:forward), fdjtype=fdtype, fdhtype=Val(:hcentral),
        relstep=nothing, absstep=nothing, dir=true
    )

# Fields

  - `fdtype::T1`: finite difference type
  - `fdjtype::T2`: finite difference type for the Jacobian
  - `fdhtype::T3`: finite difference type for the Hessian
  - `relstep`: relative finite difference step size
  - `absstep`: absolute finite difference step size
  - `dir`: direction of the finite difference step
"""
Base.@kwdef struct AutoFiniteDiff{T1, T2, T3, S1, S2, S3} <: AbstractADType
    fdtype::T1 = Val(:forward)
    fdjtype::T2 = fdtype
    fdhtype::T3 = Val(:hcentral)
    relstep::S1 = nothing
    absstep::S2 = nothing
    dir::S3 = true
end

mode(::AutoFiniteDiff) = ForwardMode()

function Base.show(io::IO, backend::AutoFiniteDiff)
    print(io, AutoFiniteDiff, "(")
    backend.fdtype != Val(:forward) &&
        print(io, "fdtype=", repr(backend.fdtype; context = io), ", ")
    backend.fdjtype != backend.fdtype &&
        print(io, "fdjtype=", repr(backend.fdjtype; context = io), ", ")
    backend.fdhtype != Val(:hcentral) &&
        print(io, "fdhtype=", repr(backend.fdhtype; context = io), ", ")
    !isnothing(backend.relstep) &&
        print(io, "relstep=", repr(backend.relstep; context = io), ", ")
    !isnothing(backend.absstep) &&
        print(io, "absstep=", repr(backend.absstep; context = io), ", ")
    backend.dir != true &&
        print(io, "dir=", repr(backend.dir; context = io))
    print(io, ")")
end

"""
    AutoFiniteDifferences{T}

Struct used to select the [FiniteDifferences.jl](https://github.com/JuliaDiff/FiniteDifferences.jl) backend for automatic differentiation.

Defined by [ADTypes.jl](https://github.com/SciML/ADTypes.jl).

# Constructors

    AutoFiniteDifferences(; fdm)

# Fields

  - `fdm::T`: a [`FiniteDifferenceMethod`](https://juliadiff.org/FiniteDifferences.jl/stable/pages/api/#FiniteDifferences.FiniteDifferenceMethod)
"""
Base.@kwdef struct AutoFiniteDifferences{T} <: AbstractADType
    fdm::T
end

mode(::AutoFiniteDifferences) = ForwardMode()

function Base.show(io::IO, backend::AutoFiniteDifferences)
    print(io, AutoFiniteDifferences, "(fdm=", repr(backend.fdm; context = io), ")")
end

"""
    AutoForwardDiff{chunksize,T}

Struct used to select the [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) backend for automatic differentiation.

Defined by [ADTypes.jl](https://github.com/SciML/ADTypes.jl).

# Constructors

    AutoForwardDiff(; chunksize=nothing, tag=nothing)

# Type parameters

  - `chunksize`: the preferred [chunk size](https://juliadiff.org/ForwardDiff.jl/stable/user/advanced/#Configuring-Chunk-Size) to evaluate several derivatives at once

# Fields

  - `tag::T`: a [custom tag](https://juliadiff.org/ForwardDiff.jl/release-0.10/user/advanced.html#Custom-tags-and-tag-checking-1) to handle nested differentiation calls (usually not necessary)
"""
struct AutoForwardDiff{chunksize, T} <: AbstractADType
    tag::T
end

AutoForwardDiff{chunksize}(tag::T) where {chunksize, T} = AutoForwardDiff{chunksize, T}(tag)

function AutoForwardDiff(; chunksize = nothing, tag = nothing)
    return AutoForwardDiff{chunksize}(tag)
end

mode(::AutoForwardDiff) = ForwardMode()

function Base.show(io::IO, backend::AutoForwardDiff{chunksize}) where {chunksize}
    print(io, AutoForwardDiff, "(")
    chunksize !== nothing && print(io, "chunksize=", repr(chunksize; context = io),
        (backend.tag !== nothing ? ", " : ""))
    backend.tag !== nothing && print(io, "tag=", repr(backend.tag; context = io))
    print(io, ")")
end

"""
    AutoTaylorDiff{order}

Struct used to select the [TaylorDiff.jl](https://github.com/JuliaDiff/TaylorDiff.jl) backend for automatic differentiation.

Defined by [ADTypes.jl](https://github.com/SciML/ADTypes.jl).

# Constructors

    AutoTaylorDiff(; order = 1)

# Type parameters

  - `order`: the order of the Taylor-mode automatic differentiation
"""
struct AutoTaylorDiff{order} <: AbstractADType end

function AutoTaylorDiff(; order = 1)
    return AutoTaylorDiff{order}()
end

mode(::AutoTaylorDiff) = ForwardMode()

function Base.show(io::IO, ::AutoTaylorDiff{order}) where {order}
    print(io, AutoTaylorDiff, "(")
    print(io, "order=", repr(order; context = io))
    print(io, ")")
end

"""
    AutoGTPSA{D}

Struct used to select the [GTPSA.jl](https://github.com/bmad-sim/GTPSA.jl) backend for automatic differentiation.

Defined by [ADTypes.jl](https://github.com/SciML/ADTypes.jl).

# Constructors

    AutoGTPSA(; descriptor=nothing)

# Fields

  - `descriptor::D`: can be either

      + a GTPSA `Descriptor` specifying the number of variables/parameters, parameter
        order, individual variable/parameter truncation orders, and maximum order. See
        the [GTPSA.jl documentation](https://bmad-sim.github.io/GTPSA.jl/stable/man/b_descriptor/) for more details.
      + `nothing` to automatically use a `Descriptor` given the context.
"""
Base.@kwdef struct AutoGTPSA{D} <: AbstractADType
    descriptor::D = nothing
end

mode(::AutoGTPSA) = ForwardMode()

function Base.show(io::IO, backend::AutoGTPSA{D}) where {D}
    print(io, AutoGTPSA, "(")
    D != Nothing && print(io, "descriptor=", repr(backend.descriptor; context = io))
    print(io, ")")
end

"""
    AutoMooncake

Struct used to select the [Mooncake.jl](https://github.com/compintell/Mooncake.jl) backend for automatic differentiation in reverse mode.

Defined by [ADTypes.jl](https://github.com/SciML/ADTypes.jl).

!!! info

    When forward mode became available in Mooncake.jl v0.4.147, another struct called [`AutoMooncakeForward`](@ref) was introduced.
    It was kept separate to avoid a breaking release of ADTypes.jl.
    [`AutoMooncake`](@ref) remains for reverse mode only.

# Constructors

    AutoMooncake(; config=nothing)

# Fields

  - `config`: either `nothing` or an instance of `Mooncake.Config` -- see the docstring of `Mooncake.Config` for more information. `AutoMooncake(; config=nothing)` is equivalent to `AutoMooncake(; config=Mooncake.Config())`, i.e. the default configuration.
"""
Base.@kwdef struct AutoMooncake{Tconfig} <: AbstractADType
    config::Tconfig = nothing
end

mode(::AutoMooncake) = ReverseMode()

function Base.show(io::IO, backend::AutoMooncake)
    print(io, AutoMooncake, "(")
    backend.config !== nothing &&
        print(io, "config=", repr(backend.config; context = io))
    print(io, ")")
end

"""
    AutoMooncakeForward

Struct used to select the [Mooncake.jl](https://github.com/compintell/Mooncake.jl) backend for automatic differentiation in forward mode.

Defined by [ADTypes.jl](https://github.com/SciML/ADTypes.jl).

!!! info

    This struct was introduced when forward mode became available in Mooncake.jl v0.4.147.
    It was kept separate from [`AutoMooncake`](@ref) to avoid a breaking release of ADTypes.jl.
    [`AutoMooncake`](@ref) remains for reverse mode only.

# Constructors

    AutoMooncakeForward(; config=nothing)

# Fields

  - `config`: either `nothing` or an instance of `Mooncake.Config` -- see the docstring of `Mooncake.Config` for more information. `AutoForwardMooncake(; config=nothing)` is equivalent to `AutoForwardMooncake(; config=Mooncake.Config())`, i.e. the default configuration.
"""
Base.@kwdef struct AutoMooncakeForward{Tconfig} <: AbstractADType
    config::Tconfig = nothing
end

mode(::AutoMooncakeForward) = ForwardMode()

function Base.show(io::IO, backend::AutoMooncakeForward)
    print(io, AutoMooncakeForward, "(")
    backend.config !== nothing &&
        print(io, "config=", repr(backend.config; context = io))
    print(io, ")")
end

"""
    AutoPolyesterForwardDiff{chunksize,T}

Struct used to select the [PolyesterForwardDiff.jl](https://github.com/JuliaDiff/PolyesterForwardDiff.jl) backend for automatic differentiation.

Defined by [ADTypes.jl](https://github.com/SciML/ADTypes.jl).

# Constructors

    AutoPolyesterForwardDiff(; chunksize=nothing, tag=nothing)

# Type parameters

  - `chunksize`: the preferred [chunk size](https://juliadiff.org/ForwardDiff.jl/stable/user/advanced/#Configuring-Chunk-Size) to evaluate several derivatives at once

# Fields

  - `tag::T`: a [custom tag](https://juliadiff.org/ForwardDiff.jl/release-0.10/user/advanced.html#Custom-tags-and-tag-checking-1) to handle nested differentiation calls (usually not necessary)
"""
struct AutoPolyesterForwardDiff{chunksize, T} <: AbstractADType
    tag::T
end

function AutoPolyesterForwardDiff{chunksize}(tag::T) where {chunksize, T}
    return AutoPolyesterForwardDiff{chunksize, T}(tag)
end

function AutoPolyesterForwardDiff(; chunksize = nothing, tag = nothing)
    return AutoPolyesterForwardDiff{chunksize}(tag)
end

mode(::AutoPolyesterForwardDiff) = ForwardMode()

function Base.show(io::IO, backend::AutoPolyesterForwardDiff{chunksize}) where {chunksize}
    print(io, AutoPolyesterForwardDiff, "(")
    chunksize !== nothing && print(io, "chunksize=", repr(chunksize; context = io),
        (backend.tag !== nothing ? ", " : ""))
    backend.tag !== nothing && print(io, "tag=", repr(backend.tag; context = io))
    print(io, ")")
end

"""
    AutoReverseDiff{compile}

Struct used to select the [ReverseDiff.jl](https://github.com/JuliaDiff/ReverseDiff.jl) backend for automatic differentiation.

Defined by [ADTypes.jl](https://github.com/SciML/ADTypes.jl).

# Constructors

    AutoReverseDiff(; compile::Union{Val, Bool} = Val(false))

# Fields

  - `compile::Union{Val, Bool}`: whether to allow pre-recording and reusing a tape (which speeds up the differentiation process).

      + If `compile=false` or `compile=Val(false)`, a new tape must be recorded at every call to the differentiation operator.
      + If `compile=true` or `compile=Val(true)`, a tape can be pre-recorded on an example input and then reused at every differentiation call.

    The boolean version of this keyword argument is taken as the type parameter.

!!! warning

    Pre-recording a tape only captures the path taken by the differentiated function _when executed on the example input_.
    If said function has value-dependent branching behavior, reusing pre-recorded tapes can lead to incorrect results.
    In such situations, you should keep the default setting `compile=Val(false)`.
    For more details, please refer to ReverseDiff's [`AbstractTape` API documentation](https://juliadiff.org/ReverseDiff.jl/dev/api/#The-AbstractTape-API).

!!! info

    Despite what its name may suggest, the `compile` setting does not prescribe whether or not the tape is compiled with [`ReverseDiff.compile`](https://juliadiff.org/ReverseDiff.jl/dev/api/#ReverseDiff.compile) after being recorded.
    This is left as a private implementation detail.
"""
struct AutoReverseDiff{C} <: AbstractADType
    compile::Bool  # this field is left for legacy reasons

    function AutoReverseDiff(; compile::Union{Val, Bool} = Val(false))
        _compile = _unwrap_val(compile)
        return new{_compile}(_compile)
    end
end

function Base.getproperty(ad::AutoReverseDiff, s::Symbol)
    if s === :compile
        Base.depwarn(
            "`ad.compile` where `ad` is `AutoReverseDiff` has been deprecated and will be removed in v2. Instead it is available as a compile-time constant as `AutoReverseDiff{true}` or `AutoReverseDiff{false}`.",
            :getproperty)
    end
    return getfield(ad, s)
end

mode(::AutoReverseDiff) = ReverseMode()

function Base.show(io::IO, ::AutoReverseDiff{compile}) where {compile}
    print(io, AutoReverseDiff, "(")
    compile && print(io, "compile=true")
    print(io, ")")
end

"""
    AutoSymbolics

Struct used to select the [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) backend for automatic differentiation.

Defined by [ADTypes.jl](https://github.com/SciML/ADTypes.jl).

# Constructors

    AutoSymbolics()
"""
struct AutoSymbolics <: AbstractADType end

mode(::AutoSymbolics) = SymbolicMode()

"""
    AutoTapir

!!! danger

    `AutoTapir` is deprecated following a package renaming, please use [`AutoMooncake`](@ref) instead.

Struct used to select the [Tapir.jl](https://github.com/compintell/Tapir.jl) backend for automatic differentiation.

Defined by [ADTypes.jl](https://github.com/SciML/ADTypes.jl).

# Constructors

    AutoTapir(; safe_mode=true)

# Fields

  - `safe_mode::Bool`: whether to run additional checks to catch errors early.
"""
struct AutoTapir <: AbstractADType
    safe_mode::Bool
end

mode(::AutoTapir) = ReverseMode()

function Base.show(io::IO, backend::AutoTapir)
    print(io, AutoTapir, "(")
    !(backend.safe_mode) && print(io, "safe_mode=false")
    print(io, ")")
end

"""
    AutoTracker

Struct used to select the [Tracker.jl](https://github.com/FluxML/Tracker.jl) backend for automatic differentiation.

Defined by [ADTypes.jl](https://github.com/SciML/ADTypes.jl).

# Constructors

    AutoTracker()
"""
struct AutoTracker <: AbstractADType end

mode(::AutoTracker) = ReverseMode()

"""
    AutoZygote

Struct used to select the [Zygote.jl](https://github.com/FluxML/Zygote.jl) backend for automatic differentiation.

Defined by [ADTypes.jl](https://github.com/SciML/ADTypes.jl).

# Constructors

    AutoZygote()
"""
struct AutoZygote <: AbstractADType end

mode(::AutoZygote) = ReverseMode()

"""
    NoAutoDiff

Struct used to select no automatic differentiation.

Defined by [ADTypes.jl](https://github.com/SciML/ADTypes.jl).

# Constructors

    NoAutoDiff()
"""
struct NoAutoDiff <: AbstractADType end

"""
    NoAutoDiffSelectedError <: Exception

Signifies that code tried to use automatic differentiation, but [`NoAutoDiff`](@ref) was specified.

# Constructor

    NoAutoDiffSelectedError(msg::String)
"""
struct NoAutoDiffSelectedError <: Exception
    msg::String
end

NoAutoDiffSelectedError() = NoAutoDiffSelectedError("Automatic differentiation can not be used with NoAutoDiff()")

function mode(::NoAutoDiff)
    throw(NoAutoDiffSelectedError())
end

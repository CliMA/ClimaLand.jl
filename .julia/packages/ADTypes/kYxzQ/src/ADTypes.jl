"""
    ADTypes.jl

[ADTypes.jl](https://github.com/SciML/ADTypes.jl) is a multi-valued logic system to choose an automatic differentiation (AD) package and specify its parameters.
"""
module ADTypes

using Base: @deprecate

"""
    AbstractADType

Abstract supertype for all AD choices.
"""
abstract type AbstractADType end

Base.broadcastable(ad::AbstractADType) = Ref(ad)

@inline _unwrap_val(::Val{T}) where {T} = T
@inline _unwrap_val(x) = x

include("compat.jl") # @public macro
include("mode.jl")
include("dense.jl")
include("sparse.jl")
include("legacy.jl")
include("symbols.jl")

# Automatic Differentiation
export AbstractADType
export AutoChainRules,
       AutoDiffractor,
       AutoEnzyme,
       AutoFastDifferentiation,
       AutoFiniteDiff,
       AutoFiniteDifferences,
       AutoForwardDiff,
       AutoGTPSA,
       AutoModelingToolkit,
       AutoMooncake,
       AutoMooncakeForward,
       AutoPolyesterForwardDiff,
       AutoReverseDiff,
       AutoSymbolics,
       AutoTapir,
       AutoTaylorDiff,
       AutoTracker,
       AutoZygote,
       NoAutoDiff,
       NoAutoDiffSelectedError
@public AbstractMode
@public ForwardMode, ReverseMode, ForwardOrReverseMode, SymbolicMode
@public mode
@public Auto

# Sparse Automatic Differentiation
export AutoSparse
@public dense_ad

# Sparsity detection
export AbstractSparsityDetector
export jacobian_sparsity, hessian_sparsity
@public sparsity_detector
@public NoSparsityDetector
@public KnownJacobianSparsityDetector
@public KnownHessianSparsityDetector

# Matrix coloring
export AbstractColoringAlgorithm
export column_coloring, row_coloring, symmetric_coloring
@public coloring_algorithm
@public NoColoringAlgorithm

# legacy exports are taken care of by @deprecated

end

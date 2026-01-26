module StatsAPI

using LinearAlgebra

include("statisticalmodel.jl")
include("regressionmodel.jl")

"""
    params(model)

Return all parameters of a model.
"""
function params end

function params! end

#    pairwise(f, x[, y])
#
# Return a matrix holding the result of applying `f` to all possible pairs
# of entries in iterators `x` and `y`, and return it.
#
# This generic function is owned by StatsBase.jl, which is the sole provider
# of the default definition.
function pairwise end

#    pairwise!(f, dest::AbstractMatrix, x[, y])
#
# Store in matrix `dest` the result of applying `f` to all possible pairs
# of entries in iterators `x` and `y`, and return it.
#
# This generic function is owned by StatsBase.jl, which is the sole provider
# of the default definition.
function pairwise! end

"""
    HypothesisTest

Abstract supertype for all statistical hypothesis tests.
Subtypes must implement [`pvalue`](@ref) at a minimum and may also
implement functions such as [`confint`](@ref), [`nobs`](@ref), and
[`dof`](@ref) as appropriate.
"""
abstract type HypothesisTest end

"""
    pvalue(test)

Compute the p-value for a given significance test.
"""
function pvalue end

end # module

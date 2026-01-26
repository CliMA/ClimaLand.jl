__precompile__(true)

module NLSolversBase

using ADTypes: AbstractADType, AutoForwardDiff, AutoFiniteDiff
import DifferentiationInterface as DI
using FiniteDiff: FiniteDiff
using ForwardDiff: ForwardDiff
import Distributed: clear!
export AbstractObjective,
       NonDifferentiable,
       OnceDifferentiable,
       TwiceDifferentiable,
       TwiceDifferentiableHV,
       value,
       value!,
       value_gradient!,
       value_jacobian!,
       gradient,
       gradient!,
       jacobian,
       jacobian!,
       hessian,
       hessian!,
       value!!,
       value_gradient!!,
       value_jacobian!!,
       hessian!!,
       hv_product,
       hv_product!,
       only_fg!,
       only_fgh!,
       only_fj!,
       only_fg,
       only_fj,
       only_g_and_fg,
       only_j_and_fj,
       only_fg_and_hv!,
       only_fghv!,
       clear!,
       f_calls,
       g_calls,
       h_calls,
       hv_calls

export AbstractConstraints, OnceDifferentiableConstraints,
    TwiceDifferentiableConstraints, ConstraintBounds

function finitediff_fdtype(autodiff)
    if autodiff == :finiteforward
        fdtype = Val{:forward}
    elseif autodiff == :finitecomplex
        fdtype = Val{:complex}
    elseif any(autodiff .== (:finite, :central, :finitecentral))
        fdtype = Val{:central}
    end
    fdtype
end

forwarddiff_chunksize(::Nothing) = nothing
forwarddiff_chunksize(::ForwardDiff.Chunk{C}) where {C} = C

is_finitediff(autodiff) = autodiff ∈ (:central, :finite, :finiteforward, :finitecomplex)
is_forwarddiff(autodiff) = autodiff ∈ (:forward, :forwarddiff, true)

get_adtype(autodiff::AbstractADType, chunk=nothing) = autodiff

function get_adtype(autodiff::Union{Symbol,Bool}, chunk=nothing)
    if is_finitediff(autodiff)
        return AutoFiniteDiff(; fdtype=finitediff_fdtype(autodiff)())
    elseif is_forwarddiff(autodiff)
        return AutoForwardDiff(; chunksize=forwarddiff_chunksize(chunk))
    else
        error("The autodiff value $autodiff is not supported. Use :finite or :forward.")
    end
end

x_of_nans(x, Tf=eltype(x)) = fill!(Tf.(x), Tf(NaN))

include("objective_types/inplace_factory.jl")
include("objective_types/abstract.jl")
include("objective_types/nondifferentiable.jl")
include("objective_types/oncedifferentiable.jl")
include("objective_types/twicedifferentiable.jl")
include("objective_types/twicedifferentiablehv.jl")
include("objective_types/incomplete.jl")
include("objective_types/constraints.jl")
include("interface.jl")

NonDifferentiable(f::OnceDifferentiable, x::AbstractArray) = NonDifferentiable(f.f, x, copy(f.F))
NonDifferentiable(f::TwiceDifferentiable, x::AbstractArray) = NonDifferentiable(f.f, x, copy(f.F))
NonDifferentiable(f::TwiceDifferentiableHV, x::AbstractArray) = NonDifferentiable(f.f, x, copy(f.F))
end # module

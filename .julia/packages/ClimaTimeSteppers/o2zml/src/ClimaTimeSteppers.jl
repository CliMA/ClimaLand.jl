"""
    ClimaTimeSteppers

Ordinary differential equation solvers

JuliaDiffEq terminology:

* _Function_: the right-hand side function df/dt.
  * by default, a function gets wrapped in an `ODEFunction`
  * define new `IncrementingODEFunction` to support incrementing function calls.

* _Problem_: Function, initial u, time span, parameters and options

  `du/dt = f(u,p,t) = fL(u,p,t)  + fR(u,p,t)`

  `fR(u,p,t) == f(u.p,t) - fL(u,p,t)`
  `fL(u,_,_) == A*u for some `A` (matrix free)`

  SplitODEProlem(fL, fR)


  * `ODEProblem` from SciMLBase.jl
    - use `jac` option to `ODEFunction` for linear + full IMEX (https://docs.sciml.ai/latest/features/performance_overloads/#ode_explicit_jac-1)
  * `SplitODEProblem` for linear + remainder IMEX
  * `MultirateODEProblem` for true multirate

* _Algorithm_: small objects (often singleton) which indicate what algorithm + options (e.g. linear solver type)
  * define new abstract `DistributedODEAlgorithm`, algorithms in this pacakge will be subtypes of this
  * define new `Multirate` for multirate solvers

* _Integrator_: contains everything necessary to solve. Used as:

  * define new `DistributedODEIntegrator` for solvers in this package

      init(prob, alg, options...) => integrator
      step!(int) => runs single step
      solve!(int) => runs it to end
      solve(prob, alg, options...) => init + solve!

* _Solution_ (not implemented): contains the "solution" to the ODE.


"""
module ClimaTimeSteppers


using KernelAbstractions
using LinearAlgebra
using LinearOperators
using StaticArrays
import ClimaComms
using NVTX

export AbstractAlgorithmName, AbstractAlgorithmConstraint, Unconstrained, SSP

array_device(::Union{Array, SArray, MArray}) = CPU()
array_device(x) = CUDADevice() # assume CUDA

import DiffEqBase, SciMLBase, LinearAlgebra, Krylov

include(joinpath("utilities", "sparse_coeffs.jl"))
include(joinpath("utilities", "fused_increment.jl"))
include("sparse_containers.jl")
include("functions.jl")

abstract type DistributedODEAlgorithm <: DiffEqBase.AbstractODEAlgorithm end

abstract type AbstractAlgorithmName end

"""
    AbstractAlgorithmConstraint

A mechanism for constraining which operations can be performed by an algorithm
for solving ODEs.

For example, an unconstrained algorithm might compute a Runge-Kutta stage by
taking linear combinations of tendencies; i.e., by adding quantities of the form
`dt * tendency(state)`. On the other hand, a "strong stability preserving"
algorithm can only take linear combinations of "incremented states"; i.e., it
only adds quantities of the form `state + dt * coefficient * tendency(state)`.
"""
abstract type AbstractAlgorithmConstraint end

"""
    Unconstrained

Indicates that an algorithm may perform any supported operations.
"""
struct Unconstrained <: AbstractAlgorithmConstraint end

default_constraint(::AbstractAlgorithmName) = Unconstrained()

"""
    SSP

Indicates that an algorithm must be "strong stability preserving", which makes
it easier to guarantee that the algorithm will preserve monotonicity properties
satisfied by the initial state. For example, this ensures that the algorithm
will be able to use limiters in a mathematically consistent way.
"""
struct SSP <: AbstractAlgorithmConstraint end

SciMLBase.allowscomplex(alg::DistributedODEAlgorithm) = true
include("integrators.jl")

include("utilities/update_signal_handler.jl")
include("utilities/convergence_condition.jl")
include("utilities/convergence_checker.jl")
include("nl_solvers/newtons_method.jl")


n_stages_ntuple(::Type{<:NTuple{Nstages}}) where {Nstages} = Nstages
n_stages_ntuple(::Type{<:SVector{Nstages}}) where {Nstages} = Nstages

# Include concrete implementations
const SPCO = SparseCoeffs

include("solvers/imex_tableaus.jl")
include("solvers/explicit_tableaus.jl")
include("solvers/imex_ark.jl")
include("solvers/imex_ssprk.jl")
include("solvers/multirate.jl")
include("solvers/lsrk.jl")
include("solvers/mis.jl")
include("solvers/wickerskamarock.jl")
include("solvers/rosenbrock.jl")

include("Callbacks.jl")

include("arbitrary_number_types.jl")

benchmark_step(integrator, device) =
    @warn "Must load CUDA, BenchmarkTools, OrderedCollections, StatsBase, PrettyTables to trigger the ClimaTimeSteppersBenchmarkToolsExt extension"

end

export LSRK54CarpenterKennedy, LSRK144NiegemannDiehlBusch, LSRKEulerMethod

"""
    LowStorageRungeKutta2N <: DistributedODEAlgorithm

A class of low-storage Runge-Kutta algorithms, which use only one additional
copy of the state vector ``u`` (often referred to as ``2N`` schemes).

The available implementations are:
 - [`LSRKEulerMethod`](@ref)
 - [`LSRK54CarpenterKennedy`](@ref)
 - [`LSRK144NiegemannDiehlBusch`](@ref)
"""
abstract type LowStorageRungeKutta2N <: DistributedODEAlgorithm end


"""
    LowStorageRungeKutta2NTableau

Storage for the tableau of a [`LowStorageRungeKutta2N`](@ref) algorithm.
"""
struct LowStorageRungeKutta2NTableau{T <: NTuple}
    "low storage RK coefficient vector A (rhs scaling)"
    A::T
    "low storage RK coefficient vector B (rhs add in scaling)"
    B::T
    "low storage RK coefficient vector C (time scaling)"
    C::T
end
n_stages(::LowStorageRungeKutta2NTableau{T}) where {T} = n_stages_ntuple(T)

struct LowStorageRungeKutta2NIncCache{T <: LowStorageRungeKutta2NTableau, A}
    tableau::T
    du::A
end

function init_cache(prob, alg::LowStorageRungeKutta2N; kwargs...)
    # @assert prob.problem_type isa DiffEqBase.IncrementingODEProblem ||
    #     prob.f isa DiffEqBase.IncrementingODEFunction
    du = zero(prob.u0)
    return LowStorageRungeKutta2NIncCache(tableau(alg, eltype(du)), du)
end

n_stages(cache::LowStorageRungeKutta2NIncCache) = n_stages(cache.tableau)

function step_u!(int, cache::LowStorageRungeKutta2NIncCache)
    (; C, A, B) = cache.tableau
    du = cache.du
    (; u, p, t, dt) = int

    for stage in 1:n_stages(cache)
        #  du .= f(u, p, t + C[stage]*dt) .+ A[stage] .* du
        stage_time = t + C[stage] * dt
        int.sol.prob.f(du, u, p, stage_time, 1, A[stage])
        u .+= (dt * B[stage]) .* du
    end
end

# for Multirate
function init_inner(prob, outercache::LowStorageRungeKutta2NIncCache, dt)
    OffsetODEFunction(prob.f.f1, zero(dt), one(dt), zero(dt), outercache.du)
end
function update_inner!(innerinteg, outercache::LowStorageRungeKutta2NIncCache, f_slow, u, p, t, dt, stage)

    (; C, A, B) = outercache.tableau
    f_offset = innerinteg.sol.prob.f
    N = n_stages(outercache)

    τ0 = t + C[stage] * dt
    τ1 = stage == N ? t + dt : t + C[stage + 1] * dt
    f_offset.α = τ0
    innerinteg.t = zero(τ0)
    DiffEqBase.add_tstop!(innerinteg, τ1 - τ0) # TODO: verify correctness

    #  du .= f(u, p, t + C[stage]*dt) .+ A[stage] .* du
    f_slow(f_offset.x, u, p, τ0, 1, A[stage])

    C0 = C[stage]
    C1 = stage == N ? one(C[stage]) : C[stage + 1]
    f_offset.γ = B[stage] / (C1 - C0)
end



"""
    LSRK54CarpenterKennedy()

The 4th-order, 5-stage [`LowStorageRungeKutta2N`])(ref) scheme from Solution
3 of [CK1994](@cite).
"""
struct LSRK54CarpenterKennedy <: LowStorageRungeKutta2N end

function tableau(::LSRK54CarpenterKennedy, RT)
    RKA = (
        RT(0),
        RT(-567301805773 // 1357537059087),
        RT(-2404267990393 // 2016746695238),
        RT(-3550918686646 // 2091501179385),
        RT(-1275806237668 // 842570457699),
    )

    RKB = (
        RT(1432997174477 // 9575080441755),
        RT(5161836677717 // 13612068292357),
        RT(1720146321549 // 2090206949498),
        RT(3134564353537 // 4481467310338),
        RT(2277821191437 // 14882151754819),
    )

    RKC = (
        RT(0),
        RT(1432997174477 // 9575080441755),
        RT(2526269341429 // 6820363962896),
        RT(2006345519317 // 3224310063776),
        RT(2802321613138 // 2924317926251),
    )

    return LowStorageRungeKutta2NTableau(RKA, RKB, RKC)
end

"""
    LSRK144NiegemannDiehlBusch()

The 4th-order, 14-stage, [`LowStorageRungeKutta2N`])(ref) scheme of
[NDB2012](@cite) with optimized stability region
"""
struct LSRK144NiegemannDiehlBusch <: LowStorageRungeKutta2N end

function tableau(::LSRK144NiegemannDiehlBusch, RT)

    RKA = (
        RT(0),
        RT(-0.7188012108672410),
        RT(-0.7785331173421570),
        RT(-0.0053282796654044),
        RT(-0.8552979934029281),
        RT(-3.9564138245774565),
        RT(-1.5780575380587385),
        RT(-2.0837094552574054),
        RT(-0.7483334182761610),
        RT(-0.7032861106563359),
        RT(0.0013917096117681),
        RT(-0.0932075369637460),
        RT(-0.9514200470875948),
        RT(-7.1151571693922548),
    )

    RKB = (
        RT(0.0367762454319673),
        RT(0.3136296607553959),
        RT(0.1531848691869027),
        RT(0.0030097086818182),
        RT(0.3326293790646110),
        RT(0.2440251405350864),
        RT(0.3718879239592277),
        RT(0.6204126221582444),
        RT(0.1524043173028741),
        RT(0.0760894927419266),
        RT(0.0077604214040978),
        RT(0.0024647284755382),
        RT(0.0780348340049386),
        RT(5.5059777270269628),
    )

    RKC = (
        RT(0),
        RT(0.0367762454319673),
        RT(0.1249685262725025),
        RT(0.2446177702277698),
        RT(0.2476149531070420),
        RT(0.2969311120382472),
        RT(0.3978149645802642),
        RT(0.5270854589440328),
        RT(0.6981269994175695),
        RT(0.8190890835352128),
        RT(0.8527059887098624),
        RT(0.8604711817462826),
        RT(0.8627060376969976),
        RT(0.8734213127600976),
    )

    LowStorageRungeKutta2NTableau(RKA, RKB, RKC)
end


"""
    LSRKEulerMethod()

An implementation of explicit Euler method using [`LowStorageRungeKutta2N`](@ref) infrastructure.
This is mainly for debugging.
"""
struct LSRKEulerMethod <: LowStorageRungeKutta2N end
tableau(::LSRKEulerMethod, RT) = LowStorageRungeKutta2NTableau((RT(0),), (RT(1),), (RT(0),))

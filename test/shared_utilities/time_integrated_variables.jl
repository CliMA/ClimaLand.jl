import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
using Test
using ClimaLand
import ClimaLand:
    AbstractExpModel,
    TimeIntegratedVariable,
    RunningMean,
    RunningSum,
    TimeIntegral,
    time_integrated_tendency!
using ClimaLand.Domains: Column

# A minimal explicit model whose prognostic variables are entirely a set of
# `TimeIntegratedVariable`s, declared through the shared-utility helpers.
struct TIVTestModel{FT, D, S} <: AbstractExpModel{FT}
    domain::D
    specs::S
end
ClimaLand.name(::TIVTestModel) = :tiv
ClimaLand.prognostic_vars(m::TIVTestModel) =
    ClimaLand.time_integrated_prognostic_vars(m.specs)
ClimaLand.prognostic_types(m::TIVTestModel) =
    ClimaLand.time_integrated_prognostic_types(m.specs)
ClimaLand.prognostic_domain_names(m::TIVTestModel) =
    ClimaLand.time_integrated_prognostic_domain_names(m.specs)

# Read the single value out of a surface field on a Column (0-D horizontal space).
scalar(field) = Array(parent(field))[1]

function build(specs; FT = Float64)
    column = Column(; zlim = (FT(-1), FT(0)), nelements = 2)
    m = TIVTestModel{FT, typeof(column), typeof(specs)}(column, specs)
    Y, p, _ = initialize(m)
    return m, Y, p
end

# Integrate `m.specs` in place with explicit Euler, returning the time series of
# the (scalar) variable named `var`. `forcing!` may mutate `p` before each step.
function euler!(m, Y, p, var, Δt, N; t0 = 0.0)
    Xf = getproperty(Y.tiv, var)
    dY = similar(Y)
    dXf = getproperty(dY.tiv, var)
    FT = eltype(parent(Xf))
    series = Vector{FT}(undef, N)
    t = FT(t0)
    for n in 1:N
        time_integrated_tendency!(dY.tiv, Y.tiv, Y, p, t, m.specs)
        @. Xf = Xf + FT(Δt) * dXf
        series[n] = scalar(Xf)
        t += FT(Δt)
    end
    return series
end

# A model whose single time-integrated variable is NamedTuple-valued, like the
# P-model acclimation capacities `AccVars`. Its prognostic type comes from the
# spec's `element_type`, exercising the structured-field declaration path.
struct NTModel{FT, D, S} <: AbstractExpModel{FT}
    domain::D
    specs::S
end
ClimaLand.name(::NTModel) = :nt
ClimaLand.prognostic_vars(m::NTModel) =
    ClimaLand.time_integrated_prognostic_vars(m.specs)
ClimaLand.prognostic_types(m::NTModel) =
    ClimaLand.time_integrated_prognostic_types(m.specs)
ClimaLand.prognostic_domain_names(m::NTModel) =
    ClimaLand.time_integrated_prognostic_domain_names(m.specs)

@testset "Declaration path (Y, not p)" begin
    FT = Float64
    specs = (
        TimeIntegratedVariable(;
            name = :precip_annual,
            reduction = RunningSum(),
            timescale = FT(365 * 86400),
            compute_instantaneous! = (dst, Y, p, t) ->
                (@. dst = p.drivers.P_liq),
        ),
        TimeIntegratedVariable(;
            name = :T_annual,
            reduction = RunningMean(),
            timescale = FT(365 * 86400),
            compute_instantaneous! = (dst, Y, p, t) ->
                (@. dst = p.drivers.T),
        ),
    )
    m, Y, p = build(specs)
    # The variables are prognostic state (live in Y), surface fields, zero-initialized.
    @test propertynames(Y.tiv) == (:precip_annual, :T_annual)
    @test ClimaLand.prognostic_types(m) == (FT, FT)
    @test ClimaLand.prognostic_domain_names(m) == (:surface, :surface)
    @test scalar(Y.tiv.precip_annual) == 0
    @test scalar(Y.tiv.T_annual) == 0
end

@testset "RunningMean: analytic + P-model EMA equivalence" begin
    FT = Float64
    τ = FT(100.0)
    f = FT(2.0)
    X0 = FT(0.5)
    spec = TimeIntegratedVariable(;
        name = :X,
        reduction = RunningMean(),
        timescale = τ,
        compute_instantaneous! = (dst, Y, p, t) -> (@. dst = p.drivers.f),
    )
    m, Y, p0 = build((spec,))
    ff = similar(Y.tiv.X)
    ff .= f
    p = (; drivers = (; f = ff))
    Y.tiv.X .= X0

    Δt = FT(1.0)
    N = 250
    series = euler!(m, Y, p, :X, Δt, N)
    X_num = series[end]

    # Explicit-Euler of dX/dt=(f-X)/τ is exactly the P-model blend X ← αX+(1-α)f.
    α = 1 - Δt / τ
    X_blend = X0
    for _ in 1:N
        X_blend = α * X_blend + (1 - α) * f
    end
    @test X_num ≈ X_blend rtol = 1e-12

    # ...which approximates the closed-form ODE solution.
    X_exact = f + (X0 - f) * exp(-N * Δt / τ)
    @test X_num ≈ X_exact rtol = 2e-2
end

@testset "RunningSum: analytic, steady state, = τ·mean" begin
    FT = Float64
    τ = FT(100.0)
    f = FT(2.0)
    spec = TimeIntegratedVariable(;
        name = :X,
        reduction = RunningSum(),
        timescale = τ,
        compute_instantaneous! = (dst, Y, p, t) -> (@. dst = p.drivers.f),
    )
    m, Y, _ = build((spec,))
    ff = similar(Y.tiv.X)
    ff .= f
    p = (; drivers = (; f = ff))
    Y.tiv.X .= 0

    Δt = FT(1.0)
    N = 250
    series = euler!(m, Y, p, :X, Δt, N)
    X_exact = τ * f * (1 - exp(-N * Δt / τ))
    @test series[end] ≈ X_exact rtol = 2e-2

    # long-run steady state → τ·f
    euler!(m, Y, p, :X, Δt, 2000)
    @test scalar(Y.tiv.X) ≈ τ * f rtol = 1e-3

    # RunningSum is exactly τ × RunningMean (same forcing, both start at 0)
    mean_spec = TimeIntegratedVariable(;
        name = :X,
        reduction = RunningMean(),
        timescale = τ,
        compute_instantaneous! = (dst, Y, p, t) -> (@. dst = p.drivers.f),
    )
    m2, Y2, _ = build((mean_spec,))
    ff2 = similar(Y2.tiv.X)
    ff2 .= f
    p2 = (; drivers = (; f = ff2))
    Y2.tiv.X .= 0
    m3, Y3, _ = build((spec,))
    ff3 = similar(Y3.tiv.X)
    ff3 .= f
    p3 = (; drivers = (; f = ff3))
    Y3.tiv.X .= 0
    mean_series = euler!(m2, Y2, p2, :X, Δt, N)
    int_series = euler!(m3, Y3, p3, :X, Δt, N)
    @test int_series[end] ≈ τ * mean_series[end] rtol = 1e-12
end

@testset "TimeIntegral: pure accumulator" begin
    FT = Float64
    f = FT(2.0)
    X0 = FT(1.0)
    spec = TimeIntegratedVariable(;
        name = :X,
        reduction = TimeIntegral(),
        timescale = FT(1.0),  # ignored
        compute_instantaneous! = (dst, Y, p, t) -> (@. dst = p.drivers.f),
    )
    m, Y, _ = build((spec,))
    ff = similar(Y.tiv.X)
    ff .= f
    p = (; drivers = (; f = ff))
    Y.tiv.X .= X0

    Δt = FT(0.5)
    N = 100
    series = euler!(m, Y, p, :X, Δt, N)
    @test series[end] ≈ X0 + f * N * Δt rtol = 1e-12  # exact for constant f
end

@testset "RunningMean is a smooth low-pass filter (sinusoidal forcing)" begin
    FT = Float64
    τ = FT(50.0)
    A = FT(3.0)
    period = FT(40.0)
    ω = FT(2π) / period
    spec = TimeIntegratedVariable(;
        name = :X,
        reduction = RunningMean(),
        timescale = τ,
        compute_instantaneous! = (dst, Y, p, t) ->
            (@. dst = A * sin(ω * t)),
    )
    m, Y, _ = build((spec,))
    p = (;)
    Y.tiv.X .= 0

    Δt = FT(0.05)
    N = round(Int, 600 / Δt)  # integrate past the transient (~5τ)
    series = euler!(m, Y, p, :X, Δt, N)

    # steady oscillation amplitude matches the first-order transfer function
    nper = round(Int, period / Δt)
    window = series[(end - nper + 1):end]
    amp_num = (maximum(window) - minimum(window)) / 2
    amp_expected = A / sqrt(1 + (ω * τ)^2)
    @test amp_num ≈ amp_expected rtol = 5e-2
    @test amp_num < A  # genuinely attenuated (smoothed)
end

@testset "A0_annual: smooth running integral vs the yearly-step scheme" begin
    FT = Float64
    year = FT(365 * 86400)
    Δt = FT(86400)              # daily, matching the current finalize cadence
    nyears = 4
    N = nyears * 365
    f0 = FT(1e-5)              # ~ potential GPP magnitude, mol C m^-2 s^-1
    trend = FT(0.15)          # +15%/yr, so annual totals differ year to year
    # synthetic daily-mean potential GPP: seasonal, non-negative, mild trend
    inst(t) = f0 * (1 + trend * (t / year)) * max(FT(0), sin(FT(2π) * t / year))

    # --- new scheme: RunningSum in Y, τ = 1 yr ---
    spec = TimeIntegratedVariable(;
        name = :A0_annual,
        reduction = RunningSum(),
        timescale = year,
        compute_instantaneous! = (dst, Y, p, t) -> (@. dst = inst(t)),
    )
    m, Y, _ = build((spec,))
    p = (;)
    Y.tiv.A0_annual .= 0
    tiv_series = euler!(m, Y, p, :A0_annual, Δt, N)

    # --- current scheme: accumulate daily, snap the 365-day sum, reset ---
    cur_series = Vector{FT}(undef, N)
    daily_totals = FT[]
    annual = FT(0)
    days_since_reset = 0
    for n in 1:N
        t = FT((n - 1) * Δt)
        push!(daily_totals, inst(t) * Δt)
        days_since_reset += 1
        if days_since_reset >= 365
            annual = sum(daily_totals)   # discrete once-a-year jump
            empty!(daily_totals)
            days_since_reset = 0
        end
        cur_series[n] = annual
    end

    tol = FT(1e-6)
    tiv_moves = count(x -> abs(x) > tol, diff(tiv_series))
    cur_moves = count(x -> abs(x) > tol, diff(cur_series))
    max_jump_tiv = maximum(abs.(diff(tiv_series)))
    max_jump_cur = maximum(abs.(diff(cur_series)))

    # the running integral updates essentially every timestep...
    @test tiv_moves >= N - 5
    # ...while the current scheme changes only ~once per year (nyears jumps)
    @test cur_moves <= nyears + 1
    # ...and those yearly jumps are far larger than any single TIV step
    @test max_jump_tiv < max_jump_cur / 10

    # both estimate the same annual total (running integral ≈ τ·mean(f) ≈ boxcar sum)
    tiv_lastyear = tiv_series[(end - 364):end]
    mean_tiv = sum(tiv_lastyear) / length(tiv_lastyear)
    @test mean_tiv ≈ cur_series[end] rtol = 0.25
    @test 50 < mean_tiv < 300   # order-of-magnitude sanity (mol C m^-2 yr^-1)
end

@testset "RunningMean on a NamedTuple-valued field (P-model AccVars analog)" begin
    FT = Float64
    col = Column(; zlim = (FT(-1), FT(0)), nelements = 2)
    fa, fb, τ = FT(5.0), FT(8.0), FT(100.0)
    accvars_type = NamedTuple{(:a, :b), Tuple{FT, FT}}
    spec = TimeIntegratedVariable(;
        name = :X,
        reduction = RunningMean(),
        timescale = τ,
        element_type = accvars_type,
        compute_instantaneous! = (dst, Y, p, t) ->
            (dst.a .= fa; dst.b .= fb; nothing),
    )
    m = NTModel{FT, typeof(col), typeof((spec,))}(col, (spec,))
    # the structured element type is carried through the declaration helper
    @test ClimaLand.prognostic_types(m) == (accvars_type,)
    Y, p, _ = initialize(m)
    Y.nt.X.a .= FT(1.0)
    Y.nt.X.b .= FT(2.0)
    dY = similar(Y)
    Δt = FT(1.0)
    N = 300
    for _ in 1:N
        time_integrated_tendency!(dY.nt, Y.nt, Y, p, FT(0), (spec,))
        Y.nt.X.a .= Y.nt.X.a .+ Δt .* dY.nt.X.a
        Y.nt.X.b .= Y.nt.X.b .+ Δt .* dY.nt.X.b
    end
    a = Array(parent(Y.nt.X.a))[1]
    b = Array(parent(Y.nt.X.b))[1]
    # each component is an independent scalar EMA toward its target
    @test a ≈ fa + (FT(1.0) - fa) * exp(-N * Δt / τ) rtol = 2e-2
    @test b ≈ fb + (FT(2.0) - fb) * exp(-N * Δt / τ) rtol = 2e-2
end

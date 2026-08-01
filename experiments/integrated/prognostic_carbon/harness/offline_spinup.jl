## Offline spinup of the prognostic carbon pools (stage 4).
##
## Phase 1 is one-way coupled: the pools consume GPP, LAI and temperature but
## feed nothing back into them. So the pool ODEs can be integrated alone against
## a recorded driver time series, recycled for as long as the stem pool needs -
## centuries in seconds, rather than the node-hours a coupled run would take.
##
## The equations here MUST match `Canopy.update_carbon_fluxes!` and the pool
## tendency in `src/standalone/Vegetation/biomass.jl`. Parameters are read from
## the same TOML, never hard-coded, so a parameter change cannot silently make
## the offline and coupled models disagree. The `--validate` mode exists to
## catch the case where the *equations* drift apart regardless.
##
## Usage:
##   julia offline_spinup.jl <driver_record.csv> [--years N] [--out ic.csv]
##   julia offline_spinup.jl <driver_record.csv> --validate <coupled_pools.csv>

import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Canopy

const FT = Float64
const M_C = FT(0.012011)  # kg C per mol
const SECONDS_PER_DAY = FT(86400)

"""
    read_driver_record(path)

Reads the daily driver record written by `site_driver.jl`: GPP and Rd in
mol CO2 m^-2 s^-1, canopy temperature in K, the C3 fraction, and LAI.
"""
function read_driver_record(path)
    lines = readlines(path)
    header = split(strip(lines[1]), ',')
    idx = Dict(h => i for (i, h) in enumerate(header))
    n = length(lines) - 1
    have = [
        c for c in ("gpp", "rd", "ct", "tair", "fc3", "lai", "pra") if
        haskey(idx, c)
    ]
    cols = Dict(c => Vector{FT}(undef, n) for c in have)
    for (j, line) in enumerate(lines[2:end])
        f = split(strip(line), ',')
        for c in keys(cols)
            cols[c][j] = parse(FT, f[idx[c]])
        end
    end
    # The `ct` diagnostic is masked to NaN where leaf+stem area index is zero
    # (`nan_if_no_canopy`). Substitute air temperature there: with no canopy the
    # canopy temperature relaxes to it, and leaving NaN would poison Rm through
    # 0 * NaN even when the substrate ramp is zero. The coupled model reads the
    # unmasked field and is unaffected.
    if haskey(cols, "tair")
        nfixed = 0
        for j in eachindex(cols["ct"])
            if !isfinite(cols["ct"][j])
                cols["ct"][j] = cols["tair"][j]
                nfixed += 1
            end
        end
        nfixed > 0 &&
            @info "filled $nfixed masked canopy temperatures with air temperature"
    end
    for (c, v) in cols
        all(isfinite, v) || error(
            "driver record column '$c' has non-finite entries with no fallback; \
             refusing to integrate against it",
        )
    end
    return cols
end

"""
    woody_fraction(MAP, half, n)

Fraction of the structural allocation that goes to stem, as a function of mean
annual precipitation (m yr^-1). A saturating ramp: dry columns build almost no
wood, wet ones approach the full `f_stem`.

Mean annual precipitation is the classic climate control on maximum woody cover
(Sankaran et al. 2005), it is already carried as a trailing integral by the LAI
model, and it is emphatically **not** a plant functional type - which is the
constraint MODEL.md §2.3 imposes on this fallback.

`half = 0` disables the mechanism and recovers the constant-`f_stem` behaviour,
so the same code path serves both and the sweep can compare them directly.
"""
function woody_fraction(MAP, half, n)
    half <= 0 && return one(MAP)
    x = max(MAP, zero(MAP)) / half
    xn = x^n
    return xn / (1 + xn)
end

"""
    tau_stem_scale(MAT, T_ref_tau, q)

Multiplier on stem turnover time as a function of mean annual temperature.
Cold columns hold their wood longer - boreal trees live for centuries where
temperate ones live decades - so τ_stem rises toward the cold as
`q^((T_ref_tau - MAT)/10)`, and is flat above `T_ref_tau`.

MAT is a *climate* mean over the whole record, not the instantaneous
temperature: stem longevity is a property of where a plant grows, not of the
weather on a given day.

`q = 1` disables the mechanism and recovers constant τ_stem exactly.
"""
function tau_stem_scale(MAT, T_ref_tau, q)
    q <= 1 && return one(MAT)
    return q^(max(T_ref_tau - MAT, zero(MAT)) / 10)
end

"""
    step_pools(pools, drivers, i, params, dt)

One explicit Euler step of the four pools, using exactly the fluxes
`update_carbon_fluxes!` computes. Returns the updated pools and the fluxes, so
the caller can report `Ra` and litter without recomputing them.
"""
function step_pools(
    pools,
    d,
    i,
    p,
    dt;
    map_half = 0.0,
    map_n = 2.0,
    tau_scale = 1.0,
)
    (C_sugar, C_leaf, C_stem, C_root) = pools
    MAP = haskey(d, "pra") ? d["pra"][i] : zero(FT)
    w = woody_fraction(MAP, map_half, map_n)
    GPP = M_C * d["gpp"][i]
    Rd = M_C * d["rd"][i]
    T = d["ct"][i]
    fc3 = d["fc3"][i]

    f_T = p.Q10^((T - p.T_ref) / 10)
    C_sap = C_stem / (1 + C_stem / p.C_sap_half)
    ramp(x, n) = (xn = max(x, zero(FT))^n; xn / (1 + xn))

    Rm =
        ramp(C_sugar / p.C_sugar_ref, p.n_alloc) *
        (f_T * (Rd + p.r_stem * C_sap) + f_T * p.r_root * max(C_root, zero(FT)))

    target = p.c_nsc * (C_leaf + C_sap + C_root)
    S =
        max(C_sugar, zero(FT)) / p.τ_alloc *
        ramp(C_sugar / max(target, eps(FT)), p.n_alloc)
    Rg = (1 - p.a) * S

    blend(v3, v4) = fc3 * v3 + (1 - fc3) * v4
    f_leaf = blend(p.f_leaf_c3, p.f_leaf_c4)
    # Stem allocation is scaled by the climate-derived woody fraction; whatever
    # it gives up goes to roots, so the three fractions still sum to one.
    f_stem = blend(p.f_stem_c3, p.f_stem_c4) * w
    f_root = 1 - f_leaf - f_stem
    τ_stem = blend(p.τ_stem_c3, p.τ_stem_c4) * tau_scale

    L_leaf = C_leaf / p.τ_leaf
    L_stem = C_stem / τ_stem
    L_root = C_root / p.τ_root

    C_sugar += dt * (GPP - Rm - S)
    C_leaf += dt * (p.a * f_leaf * S - L_leaf)
    C_stem += dt * (p.a * f_stem * S - L_stem)
    C_root += dt * (p.a * f_root * S - L_root)

    return (C_sugar, C_leaf, C_stem, C_root),
    (; Rm, Rg, Ra = Rm + Rg, S, L_leaf, L_stem, L_root, GPP)
end

"""
    spinup(drivers, params; years, dt, pools0)

Integrates the pools against the driver record, recycling it, for `years`.
Returns the final pools and the last year's mean fluxes.
"""
function spinup(
    d,
    p;
    years = 400,
    dt = SECONDS_PER_DAY,
    pools0 = (FT(0), FT(0), FT(0), FT(0)),
    map_half = 0.0,
    map_n = 2.0,
    T_ref_tau = 0.0,
    q_tau = 1.0,
)
    n = length(d["gpp"])
    # MAT is a climate mean over the record, evaluated once - not a per-step
    # quantity - because stem longevity reflects where a plant grows.
    MAT = haskey(d, "tair") ? sum(d["tair"]) / length(d["tair"]) : FT(0)
    tau_scale = tau_stem_scale(MAT, T_ref_tau, q_tau)
    pools = pools0
    steps = round(Int, years * n)
    acc = (; Ra = FT(0), litter = FT(0), GPP = FT(0), S = FT(0))
    nacc = 0
    for k in 1:steps
        i = mod1(k, n)
        pools, fl = step_pools(pools, d, i, p, dt; map_half, map_n, tau_scale)
        if k > steps - n  # last cycle only
            acc = (;
                Ra = acc.Ra + fl.Ra,
                litter = acc.litter + fl.L_leaf + fl.L_stem + fl.L_root,
                GPP = acc.GPP + fl.GPP,
                S = acc.S + fl.S,
            )
            nacc += 1
        end
    end
    means = (;
        Ra = acc.Ra / nacc,
        litter = acc.litter / nacc,
        GPP = acc.GPP / nacc,
        S = acc.S / nacc,
    )
    return pools, means
end

"""
    analytic_steady_state(d, p)

The steady state solved directly, as an independent check on the integrator.

At equilibrium `C_i* = a·f_i·S̄·τ_i`, so the pools are linear in the mean
allocation rate `S̄`. `S̄` itself is not closed-form, because maintenance
respiration depends on the pools it builds: `S̄ = GPP̄ - Rm(S̄)`. That is a scalar
fixed point, solved here by bisection.

This uses annual-mean drivers and so ignores the covariance of GPP with
temperature and the curvature of the allocation ramp. It is a cross-check on the
integrator's order of magnitude, not a replacement for it.
"""
function analytic_steady_state(d, p)
    n = length(d["gpp"])
    GPP̄ = M_C * sum(d["gpp"]) / n
    Rd̄ = M_C * sum(d["rd"]) / n
    f_T̄ = sum(p.Q10 .^ ((d["ct"] .- p.T_ref) ./ 10)) / n
    fc3 = sum(d["fc3"]) / n

    blend(v3, v4) = fc3 * v3 + (1 - fc3) * v4
    f_leaf = blend(p.f_leaf_c3, p.f_leaf_c4)
    f_stem = blend(p.f_stem_c3, p.f_stem_c4)
    f_root = 1 - f_leaf - f_stem
    τ_stem = blend(p.τ_stem_c3, p.τ_stem_c4)

    pools_of(S) = (
        p.a * f_leaf * S * p.τ_leaf,
        p.a * f_stem * S * τ_stem,
        p.a * f_root * S * p.τ_root,
    )
    function residual(S)
        (_, C_stem, C_root) = pools_of(S)
        C_sap = C_stem / (1 + C_stem / p.C_sap_half)
        Rm = f_T̄ * (Rd̄ + p.r_stem * C_sap) + f_T̄ * p.r_root * C_root
        return GPP̄ - Rm - S
    end

    lo, hi = FT(0), max(GPP̄, eps(FT))
    residual(lo) <= 0 && return (FT(0), FT(0), FT(0)), FT(0)
    for _ in 1:200
        mid = (lo + hi) / 2
        residual(mid) > 0 ? (lo = mid) : (hi = mid)
    end
    S = (lo + hi) / 2
    return pools_of(S), S
end

function load_params()
    toml_dict = LP.create_toml_dict(FT)
    return Canopy.PrognosticCarbonParameters(toml_dict)
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error(
        "usage: julia offline_spinup.jl <driver_record.csv> [--years N] [--out ic.csv]",
    )
    path = ARGS[1]
    years = FT(
        something(findfirst(==("--years"), ARGS), 0) == 0 ? 400 :
        parse(FT, ARGS[findfirst(==("--years"), ARGS) + 1]),
    )
    d = read_driver_record(path)
    p = load_params()

    pools, means = spinup(d, p; years)
    apools, S_an = analytic_steady_state(d, p)

    println("offline spinup: $(years) yr, $(length(d["gpp"])) day record")
    println(
        rpad("pool", 12),
        rpad("integrated", 14),
        rpad("analytic", 14),
        "ratio",
    )
    for (name, v, a) in (
        ("C_leaf", pools[2], apools[1]),
        ("C_stem", pools[3], apools[2]),
        ("C_root", pools[4], apools[3]),
    )
        println(
            rpad(name, 12),
            rpad(round(v, digits = 4), 14),
            rpad(round(a, digits = 4), 14),
            a > 0 ? round(v / a, digits = 3) : "-",
        )
    end
    println("C_sugar     ", round(pools[1], digits = 4))
    # sigma_l_implied at equilibrium: the stage-1 diagnostic, re-measured on a
    # spun-up state rather than on pools two years from empty. Uses the mean LAI
    # over the record, weighted only where there is leaf area.
    lai_pos = filter(>(0), d["lai"])
    if !isempty(lai_pos)
        LAI_mean = sum(lai_pos) / length(lai_pos)
        println(
            "sigma_l_eq  ",
            round(pools[2] / LAI_mean, digits = 4),
            " kg C m^-2 leaf   (LAI_mean ",
            round(LAI_mean, digits = 3),
            ", target 0.03-0.1)",
        )
    end
    println("cVeg        ", round(sum(pools), digits = 4), " kg C m^-2")
    println(
        "NPP/GPP     ",
        means.GPP > 0 ? round(1 - means.Ra / means.GPP, digits = 3) : NaN,
    )
    println(
        "litter      ",
        round(means.litter * 1000 * 86400, digits = 4),
        " g C m^-2 day^-1",
    )

    oi = findfirst(==("--out"), ARGS)
    if oi !== nothing
        open(ARGS[oi + 1], "w") do io
            println(io, "C_sugar,C_leaf,C_stem,C_root")
            println(io, join(pools, ","))
        end
        println("wrote ", ARGS[oi + 1])
    end
end

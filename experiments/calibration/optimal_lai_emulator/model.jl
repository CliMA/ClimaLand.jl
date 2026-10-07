# Offline emulator, model part. Everything is at equilibrium with the recycled 2008
# forcing: the yearly totals are the annual means of their (forcing-only) rates, and
# LAI runs hourly over the year after a 60-day warm-up.

const MM = 1 / 18.015e-3  # mol H2O m^-2 per mm

# ---- climate features (model-like definitions) ----
struct Features
    Pa::FT;
    PETa::FT;
    AI::FT;
    gd::FT;
    vpd_gs::FT;
    Tgs::FT
    A0c3a::FT;
    A0c4a::FT
    dry_all::FT    # months whose trailing 30-day P < PET/2
    dry_warm::FT   # same, counted only while the air is above freezing
    frost::FT      # mean of max(Tf - T, 0) over the year (°C)
    pcv::FT        # coefficient of variation of the monthly precipitation
end

function features(c::PointClimate)
    nh = length(c.T)
    dt = FT(3600);
    e30 = exp(-dt / (30 * DAY));
    W30 = 30 * DAY
    P30, PET30 = W30 * mean(c.precip), W30 * mean(c.PET)
    dry_all = dry_warm = 0.0
    for pass in 1:2, h in 1:nh
        P30 = P30 * e30 + W30 * c.precip[h] * (1 - e30)
        PET30 = PET30 * e30 + W30 * c.PET[h] * (1 - e30)
        if pass == 2 && P30 < PET30 / 2
            dry_all += 1
            dry_warm += c.T[h] > TF
        end
    end
    grow = c.T .> TF
    gd = 365 * mean(grow)
    vpd_gs = sum(c.VPD .* grow) / max(sum(grow), 1)
    Tgs = sum(grow) > 0 ? sum((c.T .- TF) .* grow) / sum(grow) : 0.0
    Pm = [mean(c.precip[c.month .== m]) for m in 1:12]
    Pa = YEAR * mean(c.precip);
    PETa = YEAR * mean(c.PET)
    return Features(
        Pa,
        PETa,
        PETa / max(Pa, eps()),
        gd,
        vpd_gs,
        Tgs,
        YEAR * mean(c.A0c3),
        YEAR * mean(c.A0c4),
        12 * dry_all / nh,
        12 * dry_warm / nh,
        mean(max.(TF .- c.T, 0)),
        std(Pm) / max(mean(Pm), eps()),
    )
end

# Ratio of the yearly evaporation of a bucket (capacity Wcap mm, demand PET) fed by the
# actual precipitation to that of the same bucket fed by its annual mean: below 1 where
# the rain comes when the demand is low (or in bursts that overflow the bucket).
function phase_factor(c::PointClimate, Wcap_mm)
    Wcap = Wcap_mm * MM;
    dt = FT(3600);
    P̄ = mean(c.precip)
    aet(uniform) = begin
        W = Wcap / 2;
        E = 0.0
        for pass in 1:3, h in eachindex(c.precip)
            W = min(W + dt * (uniform ? P̄ : c.precip[h]), Wcap)
            e = min(dt * c.PET[h] * W / Wcap, W)
            W -= e
            pass == 3 && (E += e)
        end
        E
    end
    return aet(false) / max(aet(true), eps())
end

# ---- configuration ----
Base.@kwdef struct Config
    z_tree::FT = LPAR.z_tree
    z_grass::FT = LPAR.z_grass
    sigma::FT = LPAR.sigma
    alpha::FT = LPAR.alpha
    f0_max::FT = LPAR.f0_max
    ai_peak::FT = 1.9
    w_humid::FT = 0.604   # width of f0(AI) on the humid side
    w_arid::FT = 0.604
    phase_exp::FT = 0.0   # f0 *= phase_factor^phase_exp
    ever_floor::FT = 0.0  # evergreen canopy keeps ever_floor × LAI_max
    m_maint::FT = 0.0     # maintenance share of the leaf cost, scaled by the season
    q10::FT = 2.0         # temperature sensitivity of that maintenance
end

# Yearly maintenance load of a canopy: days above freezing weighted by
# q10^((T - 25 °C)/10), over the year (1 for a year-round 25 °C season)
maintenance_load(c::PointClimate, q10) =
    mean(T > TF ? q10^((T - TF - 25) / 10) : 0.0 for T in c.T)
# leaf cost scale (1 - m) + m M
cost_scale(c::PointClimate, cfg) =
    cfg.m_maint == 0 ? 1.0 :
    (1 - cfg.m_maint) + cfg.m_maint * maintenance_load(c, cfg.q10)

f0_curve(AI, cfg) =
    cfg.f0_max * exp(
        -(AI < cfg.ai_peak ? cfg.w_humid : cfg.w_arid) *
        log(AI / cfg.ai_peak)^2,
    )

# LAI_max of a C3 tree canopy, with χ at the growing-season temperature and VPD
function tree_lai_max(
    c,
    f,
    cfg;
    f0 = f0_curve(f.AI, cfg),
    zs = cost_scale(c, cfg),
)
    χ = Canopy.optimal_chi(TF + f.Tgs, c.P, C_CO2, f.vpd_gs, LAIPM.β_c3, PC)
    return Canopy.compute_L_max(
        f.A0c3a,
        LPAR.k,
        cfg.z_tree * zs,
        f.Pa,
        f0,
        C_CO2 * c.P,
        χ,
        f.vpd_gs,
    )
end

# Monthly LAI over the year (and the tree share, evergreen share used)
function emulate(
    c::PointClimate,
    f::Features,
    tree::FT,
    ever::FT,
    φ::FT,
    cfg::Config,
)
    k = LPAR.k;
    nh = length(c.T);
    dt = FT(3600)
    e3 = exp(-dt / (3 * DAY));
    elai = exp(-dt * cfg.alpha / DAY)
    z = Canopy.leaf_cost(tree, cfg.z_tree, cfg.z_grass) * cost_scale(c, cfg)
    comp = Canopy.canopy_composition(
        tree,
        Canopy.open_canopy_c4_share(f.A0c3a, f.A0c4a, LPAR),
    )
    fc3 = 1 - comp.c4_grass
    A0a = fc3 * f.A0c3a + (1 - fc3) * f.A0c4a
    f0 = f0_curve(f.AI, cfg) * φ^cfg.phase_exp
    ca = C_CO2 * c.P
    χgs =
        fc3 *
        Canopy.optimal_chi(TF + f.Tgs, c.P, C_CO2, f.vpd_gs, LAIPM.β_c3, PC) +
        (1 - fc3) *
        Canopy.optimal_chi(TF + f.Tgs, c.P, C_CO2, f.vpd_gs, LAIPM.β_c4, PC)
    Lfloor =
        cfg.ever_floor * tree * ever > 0 ?
        cfg.ever_floor *
        Canopy.compute_L_max(A0a, k, z, f.Pa, f0, ca, χgs, f.vpd_gs) : 0.0
    E = tree * ever
    A0d = A0a / 365;
    L = NaN
    acc = zeros(12);
    cnt = zeros(Int, 12)
    for (pass, hs) in ((1, (nh - 60 * 24 + 1):nh), (2, 1:nh)), h in hs
        A0 = fc3 * c.A0c3[h] + (1 - fc3) * c.A0c4[h]
        χ = fc3 * c.χc3[h] + (1 - fc3) * c.χc4[h]
        A0d = A0d * e3 + DAY * A0 * (1 - e3)
        Lopt = Canopy.compute_L_steady_target(
            A0d,
            k,
            A0a,
            z,
            f.gd,
            cfg.sigma,
            f.Pa,
            f0,
            ca,
            χ,
            f.vpd_gs,
        )
        isnan(L) && (L = Lopt)
        L = L * elai + Lopt * (1 - elai)
        if pass == 2
            Lout = L + E * max(Lfloor - L, 0.0)
            m = c.month[h];
            acc[m] += Lout;
            cnt[m] += 1
        end
    end
    return acc ./ cnt
end

# ---- statistics ----
logistic(x) = 1 / (1 + exp(-x))
# weighted logistic regression of a fraction y ∈ [0, 1] (quasi-binomial, IRLS)
function fit_logit(X, y, w; iters = 50, ridge = 1e-6)
    β = zeros(size(X, 2))
    for _ in 1:iters
        μ = logistic.(X * β);
        v = max.(μ .* (1 .- μ), 1e-6)
        zz = X * β .+ (y .- μ) ./ v;
        W = w .* v
        β = (X' * (W .* X) + ridge * I) \ (X' * (W .* zz))
    end
    return β
end

wmean(x, w) = sum(x .* w) / sum(w)
function score(sims, idx; obs = OBS, w = WGT)
    S = reduce(vcat, (s' for s in sims));
    O = obs[idx, :];
    ww = w[idx]
    ann = vec(mean(S; dims = 2)) .- vec(mean(O; dims = 2))
    rmseA = sqrt(wmean(ann .^ 2, ww));
    bias = wmean(ann, ww)
    rmseM = sqrt(wmean(vec(mean((S .- O) .^ 2; dims = 2)), ww))
    nh = [PTS[i][2] > 30 for i in idx]
    djf(X) = wmean(vec(mean(X[nh, [12, 1, 2]]; dims = 2)), ww[nh])
    jja(X) = wmean(vec(mean(X[nh, [6, 7, 8]]; dims = 2)), ww[nh])
    return (;
        rmseA,
        rmseM,
        bias,
        djf = (djf(S), djf(O)),
        jja = (jja(S), jja(O)),
    )
end
fmt(r) = @sprintf(
    "LAI RMSE ann %.3f mon %.3f bias %+.3f | NH>30 DJF %.2f (obs %.2f) JJA %.2f (obs %.2f)",
    r.rmseA,
    r.rmseM,
    r.bias,
    r.djf[1],
    r.djf[2],
    r.jja[1],
    r.jja[2]
)
function share_stats(x, y, w)
    ok = isfinite.(x) .& isfinite.(y)
    x, y, w = x[ok], y[ok], w[ok]
    mx, my = wmean(x, w), wmean(y, w)
    r =
        wmean((x .- mx) .* (y .- my), w) /
        sqrt(wmean((x .- mx) .^ 2, w) * wmean((y .- my) .^ 2, w))
    return (; rmse = sqrt(wmean((x .- y) .^ 2, w)), bias = mx - my, r)
end

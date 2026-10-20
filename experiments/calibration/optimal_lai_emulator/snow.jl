# Water input of a degree-day snow store: snowfall (precipitation below freezing)
# accumulates and melts at `ddf_mm` mm per °C per day above freezing.
function water_input(c::PointClimate, ddf_mm)
    dt = FT(3600);
    ddf = ddf_mm * MM / DAY  # mol m^-2 s^-1 K^-1
    S = 0.0;
    out = similar(c.precip)
    for pass in 1:3, h in eachindex(c.T)
        if c.T[h] <= TF
            S += c.precip[h] * dt;
            w = 0.0
        else
            melt = min(S, ddf * (c.T[h] - TF) * dt);
            S -= melt
            w = c.precip[h] + melt / dt
        end
        pass == 3 && (out[h] = w)
    end
    return out
end
# moist mask (warm, and the 30-day water input ≥ PET/2) and dry warm months
function moist_mask(c::PointClimate, input)
    dt = FT(3600);
    e30 = exp(-dt / (30 * DAY));
    W30 = 30 * DAY
    P30, PET30 = W30 * mean(input), W30 * mean(c.PET);
    m = falses(length(c.T))
    for pass in 1:2, h in eachindex(c.T)
        P30 = P30 * e30 + W30 * input[h] * (1 - e30);
        PET30 = PET30 * e30 + W30 * c.PET[h] * (1 - e30)
        pass == 2 && (m[h] = c.T[h] > TF && P30 >= PET30 / 2)
    end
    return m
end
# the climate (χ at the moist-season VPD) and features with a given water input
function moist_climate(i, input)
    c, f = CLIM[i], FEAT[i]
    wet = moist_mask(c, input);
    warm = c.T .> TF
    md = 365 * mean(wet);
    w = min(md / 30, 1.0)
    vm = sum(wet) > 0 ? sum(c.VPD .* wet) / sum(wet) : 0.0
    vpd = w * vm + (1 - w) * f.vpd_gs
    dry_warm = 12 * mean(warm .& .!wet)
    c2, f2 = with_vpd(c, f, vpd)
    f3 = Features(
        f2.Pa,
        f2.PETa,
        f2.AI,
        f2.gd,
        f2.vpd_gs,
        f2.Tgs,
        f2.A0c3a,
        f2.A0c4a,
        f2.dry_all,
        dry_warm,
        f2.frost,
        f2.pcv,
    )
    return c2, f3
end

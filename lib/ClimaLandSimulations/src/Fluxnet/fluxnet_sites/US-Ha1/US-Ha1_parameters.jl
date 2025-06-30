export harvard_default_params,
    hetero_resp_harvard,
    auto_resp_harvard,
    soil_harvard,
    radiative_transfer_harvard,
    canopy_energy_balance_harvard,
    conductance_harvard,
    photosynthesis_harvard,
    plant_hydraulics_harvard

function harvard_default_params(;
    hetero_resp = hetero_resp_harvard(),
    auto_resp = auto_resp_harvard(),
    soil = soil_harvard(),
    radiative_transfer = radiative_transfer_harvard(),
    canopy_energy_balance = canopy_energy_balance_harvard(),
    conductance = conductance_harvard(),
    photosynthesis = photosynthesis_harvard(),
    plant_hydraulics = plant_hydraulics_harvard(),
)
    return SiteParams(
        hetero_resp,
        auto_resp,
        soil,
        radiative_transfer,
        canopy_energy_balance,
        conductance,
        photosynthesis,
        plant_hydraulics,
    )
end

function hetero_resp_harvard(;
    D_ref = FT(1.39e-5),
    D_liq = FT(3.17),
    α_sx = FT(194e3),
    Ea_sx = FT(61e3),
    kM_sx = FT(5e-3),
    kM_o2 = FT(0.004),
    O2_a = FT(0.209),
    D_oa = FT(1.67),
    p_sx = FT(0.024),
)
    return HeteroRespP(
        D_ref,
        D_liq,
        α_sx,
        Ea_sx,
        kM_sx,
        kM_o2,
        O2_a,
        D_oa,
        p_sx,
    )
end

function auto_resp_harvard(;
    ne = FT(8 * 1e-4),
    ηsl = FT(0.01),
    σl = FT(0.05),
    μr = FT(1.0),
    μs = FT(0.1),
    Rel = FT(0.25),
)
    return AutotrophicRespirationParameters(ne, ηsl, σl, μr, μs, Rel)
end

function soil_harvard(;
    vg_α = FT(0.04),
    vg_n = FT(2.05),
    ν = FT(0.5), # m3/m3
    ν_ss_quartz = FT(0.1),
    ν_ss_om = FT(0.1),
    K_sat = FT(4e-7), # m/s, matches Natan
    S_s = FT(1e-3), # 1/m, guess
    θ_r = FT(0.067), # m3/m3, from Wang et al. 2021 https://doi.org/10.5194/gmd-14-6741-2021
    ν_ss_gravel = FT(0.0),
    z_0m = FT(0.01),
    z_0b = FT(0.001),
    emissivity = FT(0.98),
    PAR_albedo = FT(0.2),
    NIR_albedo = FT(0.2),
)
    hydrology_cm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
    return EnergyHydrologyParameters(
        FT;
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm,
        K_sat,
        S_s,
        θ_r,
        albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(; PAR_albedo, NIR_albedo),
        emissivity,
        z_0m,
        z_0b,
    )
end

function radiative_transfer_harvard(;
    Ω = FT(0.69),
    ld = ConstantGFunction(FT(0.5)),
    α_PAR_leaf = FT(0.1),
    λ_γ_PAR = FT(5e-7),
    τ_PAR_leaf = FT(0.05),
    α_NIR_leaf = FT(0.45),
    τ_NIR_leaf = FT(0.25),
    n_layers = UInt64(20),
    ϵ_canopy = FT(0.97),
)
    return TwoStreamParameters(
        G_Function = ld,
        α_PAR_leaf = α_PAR_leaf,
        τ_PAR_leaf = τ_PAR_leaf,
        α_NIR_leaf = α_NIR_leaf,
        τ_NIR_leaf = τ_NIR_leaf,
        ϵ_canopy = ϵ_canopy,
        Ω = Ω,
        λ_γ_PAR = λ_γ_PAR,
        n_layers = n_layers,
    )
end

function canopy_energy_balance_harvard(; ac_canopy = FT(2.5e3))
    return BigLeafEnergyParameters(ac_canopy)
end

function conductance_harvard(;
    g1 = FT(141), # Wang et al: 141 sqrt(Pa) for Medlyn model; Natan used 300.
    Drel = FT(1.6),
    g0 = FT(1e-4),
)
    return MedlynConductanceParameters(Drel, g0, g1)
end

function photosynthesis_harvard(;
    is_c3 = FT(1),
    oi = FT(0.209),
    ϕ = FT(0.6),
    θj = FT(0.9),
    fC3 = FT(0.015),
    sc = FT(2e-6), # Bonan's book: range of 2-5e-6
    pc = FT(-2e6), # Bonan's book: -2e6
    Vcmax25 = FT(9e-5), # from Yujie's paper 4.5e-5 , Natan used 9e-5
    Γstar25 = FT(4.275e-5),
    Kc25 = FT(4.049e-4),
    Ko25 = FT(0.2874),
    To = FT(298.15),
    ΔHkc = FT(79430),
    ΔHko = FT(36380),
    ΔHVcmax = FT(58520),
    ΔHΓstar = FT(37830),
    ΔHJmax = FT(43540),
    ΔHRd = FT(46390),
    fC4 = FT(0.025),
    Q10 = FT(2),
    s1 = FT(0.3),
    s2 = FT(313.15),
    s3 = FT(0.2),
    s4 = FT(288.15),
    s5 = FT(1.3),
    s6 = FT(328.15),
    E = FT(0.05),
)
    return FarquharParameters(
        Vcmax25,
        Γstar25,
        Kc25,
        Ko25,
        ΔHkc,
        ΔHko,
        ΔHVcmax,
        ΔHΓstar,
        ΔHJmax,
        ΔHRd,
        To,
        oi,
        ϕ,
        θj,
        fC3,
        fC4,
        sc,
        pc,
        Q10,
        s1,
        s2,
        s3,
        s4,
        s5,
        s6,
        E,
        is_c3,
    )
end

function plant_hydraulics_harvard(;
    SAI = FT(1.0), # m2/m2 or: estimated from Wang et al, FT(0.00242) ?
    f_root_to_shoot = FT(3.5),
    K_sat_plant = 5e-9, # m/s # seems much too small?
    ψ63 = FT(-4 / 0.0098), # / MPa to m, Holtzman's original parameter value is -4 MPa
    Weibull_param = FT(4), # unitless, Holtzman's original c param value
    a = FT(0.05 * 0.0098), # Holtzman's original parameter for the bulk modulus of elasticity
    capacity = FT(10), # kg/m^2
    S_s = FT(1e-2 * 0.0098), # m3/m3/MPa to m3/m3/m
    rooting_depth = FT(0.5), # from Natan
)
    conductivity_model =
        PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
    retention_model = PlantHydraulics.LinearRetentionCurve{FT}(a)

    return PlantHydraulicsP(
        SAI,
        f_root_to_shoot,
        capacity,
        S_s,
        rooting_depth,
        conductivity_model,
        retention_model,
    )
end

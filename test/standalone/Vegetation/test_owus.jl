using Test
using ClimaLand                     # loads the package (and Canopy)
using ClimaLand.Canopy.OWUSStomata  # <-- bring the submodule into scope

@testset "OWUS smoke" begin
    owus = OWUSStomatalModel(fww=0.6, s_star=0.4, s_w=0.15)

    s_vals = 0.0:0.05:1.0
    E0     = 3.5e-3 / 86400  # m/s
    VPD    = 1500.0          # Pa
    P_air  = 101325.0        # Pa
    T_air  = 298.15          # K

    for s in s_vals
        gref = Ref(0.0)
        stomatal_conductance!(gref, owus; s=s, E0=E0, VPD=VPD, P_air=P_air, T_air=T_air)
        @test gref[] ≥ 0
    end
end


# more elaborate
# --------- helpers ----------
# expected g_sw given OWUS piecewise β(s) and the E_mol ↔ gsw relation
function expected_gsw(m::OWUSStomatalModel, s, E0, VPD, P_air)
    # piecewise β(s)
    β = if s <= m.s_w
        0.0
    elseif s <= m.s_star
        m.fww * (s - m.s_w) / max(m.s_star - m.s_w, eps())
    else
        m.fww
    end
    # T = β * E0  [m s^-1]
    ρw = 1000.0       # kg m^-3
    Mw = 0.01801528   # kg mol^-1
    E_mol = (ρw * (β * E0)) / Mw           # mol m^-2 s^-1
    gsw   = VPD > 0 ? E_mol * (P_air / VPD) : 0.0
    return min(gsw, m.gsw_max)
end

# a tiny wrapper to call the in-package function under test
function gsw_call(m::OWUSStomatalModel; s, E0, VPD, P_air, T_air=298.15)
    gref = Ref(0.0)
    stomatal_conductance!(gref, m; s=s, E0=E0, VPD=VPD, P_air=P_air, T_air=T_air)
    return gref[]
end

# canonical forcing used across tests
const E0 = 3.5e-3 / 86400        # m s^-1  (~3.5 mm/day)
const VPD = 1500.0               # Pa
const P_air = 101325.0           # Pa
const T_air = 298.15             # K
const FT = Float64

@testset "OWUS elaborate" begin
    # a model with easy-to-check shape
    m = OWUSStomatalModel(fww=0.6, s_star=0.4, s_w=0.15, gsw_max=Inf)

    @testset "Non-negativity over [0,1]" begin
        for s in 0.0:0.02:1.0
            g = gsw_call(m; s=s, E0=E0, VPD=VPD, P_air=P_air, T_air=T_air)
            @test g ≥ 0
            @test isfinite(g)
        end
    end

    @testset "Zero below s_w" begin
        for s in (0.0, 0.05, 0.149, m.s_w, m.s_w - 1e-12)
            g = gsw_call(m; s=s, E0=E0, VPD=VPD, P_air=P_air)
            @test g ≈ 0 atol=1e-14
        end
    end

    @testset "Plateau at fww for s ≥ s_star" begin
        for s in (m.s_star, 0.6, 0.8, 1.0)
            g        = gsw_call(m; s=s, E0=E0, VPD=VPD, P_air=P_air)
            g_expect = expected_gsw(m, s, E0, VPD, P_air)
            @test g ≈ g_expect rtol=1e-12
        end
    end

    @testset "Continuity at s_star" begin
        s_left  = prevfloat(m.s_star)
        s_right = nextfloat(m.s_star)
        gL = gsw_call(m; s=s_left,  E0=E0, VPD=VPD, P_air=P_air)
        gR = gsw_call(m; s=s_right, E0=E0, VPD=VPD, P_air=P_air)
        @test isapprox(gL, gR; rtol=1e-10, atol=1e-12)
    end

    @testset "Linear segment proportionality on (s_w, s_star)" begin
        s1 = 0.20   # in (0.15, 0.40)
        s2 = 0.30
        g1 = gsw_call(m; s=s1, E0=E0, VPD=VPD, P_air=P_air)
        g2 = gsw_call(m; s=s2, E0=E0, VPD=VPD, P_air=P_air)
        ratio_g = g1 / g2
        ratio_s = (s1 - m.s_w) / (s2 - m.s_w)
        @test isapprox(ratio_g, ratio_s; rtol=1e-12)
    end

    @testset "gsw_max capping" begin
        # make a huge E0 to trigger the cap
        mcap = OWUSStomatalModel(fww=0.9, s_star=0.3, s_w=0.1, gsw_max=FT(0.2)) # mol m^-2 s^-1
        s = 0.8
        g_uncapped = expected_gsw(OWUSStomatalModel(fww=0.9, s_star=0.3, s_w=0.1, gsw_max=Inf),
                                  s, E0*1e3, VPD, P_air)
        g_capped = gsw_call(mcap; s=s, E0=E0*1e3, VPD=VPD, P_air=P_air)
        @test g_uncapped > mcap.gsw_max
        @test g_capped   ≈ mcap.gsw_max rtol=0 atol=0
    end

    @testset "Builder: owus_shape_from_climaland sanity" begin
        pars = owus_shape_from_climaland(;
            E0 = 3.0e-3,        # m day^-1
            kx_max = 6e-5,      # kg m^-1 MPa^-1 s^-1 (leaf-specific)
            LAI = 2.0,
            hc = 10.0,
            Td = 86400.0,
            ks_sat = 0.5,       # m day^-1
            RAI = 1.0,
            dr = 3e-4,
            Zr = 1.0,
            psi_g50 = -2.0,     # MPa
            psi_x50 = -3.0,     # MPa
            psi_s_sat = -0.01,  # MPa
            b = 4.5
        )
        @test 0.0 ≤ pars.fww ≤ 1.0
        @test 0.0 ≤ pars.s_w  ≤ pars.s_star ≤ 1.0
    end

    @testset "Builder: build_owus_from_traits" begin
        m2 = build_owus_from_traits(;
            E0 = 3.0e-3, kx_max = 6e-5, LAI = 2.0, hc = 10.0, Td = 86400.0,
            ks_sat = 0.5, RAI = 1.0, dr = 3e-4, Zr = 1.0,
            psi_g50 = -2.0, psi_x50 = -3.0, psi_s_sat = -0.01, b = 4.5
        )
        @test m2 isa OWUSStomatalModel
        @test 0.0 ≤ m2.fww ≤ 1.0
    end

    @testset "Builder: build_owus_from_ClimaLand (mock structs)" begin
        # minimal NamedTuples that look like CL parameter structs
        canopy_params = (; LAI=2.0, canopy_height=10.0, kx_max=6e-5, psi_g50=-2.0, psi_x50=-3.0)
        soil_params   = (; Ksat=5e-6, psi_sat=-0.01, b=4.5)  # Ksat in m s^-1 (conversion inside)
        root_params   = (; RAI=1.0, dr=3e-4, Zr=1.0)
        m3 = build_owus_from_ClimaLand(;
            canopy_params=canopy_params, soil_params=soil_params, root_params=root_params
        )
        @test m3 isa OWUSStomatalModel
        # spot-check: below wilting → 0; above plateau → equals plateau
        g_lo = gsw_call(m3; s=m3.s_w - 1e-6, E0=E0, VPD=VPD, P_air=P_air)
        g_hi = gsw_call(m3; s=1.0,           E0=E0, VPD=VPD, P_air=P_air)
        @test g_lo ≈ 0 atol=1e-14
        g_plateau = expected_gsw(m3, 1.0, E0, VPD, P_air)
        @test g_hi ≈ g_plateau rtol=1e-12
    end

    @testset "Type stability (scalar call)" begin
        mT = OWUSStomatalModel(fww=FT(0.6), s_star=FT(0.4), s_w=FT(0.15))
        g = @inferred gsw_call(mT; s=FT(0.5), E0=FT(E0), VPD=FT(VPD), P_air=FT(P_air))
        @test g isa FT
    end
end

@testset "Isohydry index diagnostic" begin
    owus = OWUSStomatalModel(fww=0.6, s_star=0.4, s_w=0.15)

    # Fake time series (replace with model outputs/forcing later)
    N      = 100
    s      = range(1.0, 0.2; length=N)          # drying trajectory
    E0     = fill(2.5e-3/86400, N)              # m s^-1 (≈2.5 mm/day)
    VPD    = fill(1500.0, N)                    # Pa
    P_air  = fill(101325.0, N)                  # Pa
    T_air  = fill(298.15, N)                    # K
    ψ_soil = range(-0.2, -2.5; length=N)        # MPa

    # --- Linear closure test
    res_lin = isohydry_index(owus;
        s, E0, VPD, P_air, T_air, ψ_soil,
        hydraulics=:linear, K_lin=0.02
    )
    @test isfinite(res_lin.slope)
    @test res_lin.R2 ≥ 0 && res_lin.R2 ≤ 1

    # --- Weibull closure test
    res_weib = isohydry_index(owus;
        s, E0, VPD, P_air, T_air, ψ_soil,
        hydraulics=:weibull, Kmax=0.03, P50=-2.5, a=2.0
    )
    @test isfinite(res_weib.slope)
    @test res_weib.R2 ≥ 0 && res_weib.R2 ≤ 1

    # Optional: sanity check expected slope range
    @test res_lin.slope ≥ -0.1 && res_lin.slope ≤ 1.1
    @test res_weib.slope ≥ -0.1 && res_weib.slope ≤ 1.1
end


###############################
# scripts/calibrate_hier_eki.jl
###############################
using ClimaLand
using LinearAlgebra
using Random


# -------------------------------------------
# Load data (drivers, obs LE, site covariates)
# -------------------------------------------
# TODO: Replace with your real data loading
S = 40 # number of sites
driver_data = ... # Vector of site driver structs (each holds VPD, SWC, Rn, T, etc.)
obs_LE = ... # Vector of observed LE stacked across sites & times
site_covars = ... # (p+1) x S matrix, standardized covariates (col = site)


# -------------------------------------------
# Define hierarchy parameterization for EKI
# -------------------------------------------
K = 3 # traits (fww, s*, s_w)
p1 = size(site_covars,1) # intercept+covariates


# Flattened parameter vector: θ = vec(Γ) ⧺ vec(U)
function pack_params(Γ, U)
return vcat(vec(Γ), vec(U))
end


function unpack_params(θ, K, p1, S)
i1 = K*p1
Γ = reshape(θ[1:i1], K, p1)
U = reshape(θ[i1+1:end], K, S)
return Γ, U
end


# -------------------------------------------
# Forward model for EKI
# -------------------------------------------
function forward_LE(θ, cfg)
Γ, U = unpack_params(θ, cfg.K, cfg.p1, cfg.S)
# Precompute traits per site
traits = Vector{NTuple{3,Float64}}(undef, cfg.S)
for s in 1:cfg.S
η = Γ * cfg.w[:,s] + U[:,s]
α,a,b = η
fww = 1/(1+exp(-α))
sw = 1/(1+exp(-a))
sb = 1/(1+exp(-b))
sstar = sw + (1-sw)*sb
traits[s] = (fww,sstar,sw)
end
# Loop over obs rows
yhat = similar(cfg.obs_LE)
# U_est are site deviations; report if needed



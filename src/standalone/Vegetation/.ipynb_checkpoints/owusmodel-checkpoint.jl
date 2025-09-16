module OWUSStomata

export OWUSStomatalModel,
       stomatal_conductance!,
       owus_shape_from_climaland,
       build_owus_from_traits,
       build_owus_from_ClimaLand

"""
    OWUSStomatalModel{FT}

Stomatal model following Bassiouni et al. (2023), parameterized by
the piece-wise linear β(s) ≈ fww for s ≥ s*, linear decline to 0 at s_w.
Adapts Bassiouni et al. 2023 (https://github.com/maoyab/OWUS/tree/main) β(s) 
to a prognostic conductance usable in land-surface time stepping.

Fields
- fww::FT   : well-watered transpiration fraction (T/E₀ at high s)
- s_star::FT: soil saturation threshold where down-regulation begins (s*)
- s_w::FT   : soil saturation where transpiration ceases (s_w)
- gsw_max::FT (optional): cap on stomatal conductance (mol m⁻² s⁻¹); use Inf if unused
"""
struct OWUSStomatalModel{FT}
    fww::FT
    s_star::FT
    s_w::FT
    gsw_max::FT
end

OWUSStomatalModel(; fww::Real=0.6, s_star::Real=0.35, s_w::Real=0.10, gsw_max::Real=Inf) = begin
    FT = promote_type(typeof(fww), typeof(s_star), typeof(s_w), typeof(gsw_max), Float64)
    OWUSStomatalModel{FT}(FT(fww), FT(s_star), FT(s_w), FT(gsw_max))
end

# Physical constants (typed helpers)
const _ρw64 = 1000.0             # kg m^-3
const _Mw64 = 0.01801528         # kg mol^-1
const _R64  = 8.314462618        # J mol^-1 K^-1

@inline _consts(::Type{FT}) where {FT} = (FT(_ρw64), FT(_Mw64), FT(_R64))

# β(s) from paper (Eq. 3): 0, linear, plateau at fww
@inline function _beta_piecewise(s::FT, m::OWUSStomatalModel{FT}) where {FT}
    fww, sstar, sw = m.fww, m.s_star, m.s_w
    if s <= sw
        return zero(FT)
    elseif s <= sstar
        return fww * (s - sw) / max(sstar - sw, eps(FT))
    else
        return fww
    end
end

"""
    stomatal_conductance!(gsw_out, model; s, E0, VPD, P_air, T_air, leaf_factor=1)

Given soil saturation `s` and potential evaporation `E0` [m s^-1] at canopy,
computes the stomatal conductance (mol m^-2 s^-1) that realizes
T = β(s)*E0, via E_mol = (ρw*T)/Mw and E_mol = g_sw * (VPD/P_air).
"""
function stomatal_conductance!(
    gsw_out::Base.RefValue{FT},
    model::OWUSStomatalModel{FT};
    s::Real, E0::Real, VPD::Real, P_air::Real, T_air::Real, leaf_factor::Real=1
) where {FT}

    sFT, E0FT, VPDFT, P_airFT, T_airFT, leafFT =
        FT(s), FT(E0), FT(VPD), FT(P_air), FT(T_air), FT(leaf_factor)
    ρw, Mw, _ = _consts(FT)

    β = _beta_piecewise(sFT, model)
    T_mps = β * E0FT * leafFT                  # [m s^-1]
    E_mol = (ρw * T_mps) / Mw                  # [mol m^-2 s^-1]
    gsw   = VPDFT > eps(FT) ? E_mol * (P_airFT / VPDFT) : zero(FT)
    gsw_out[] = min(gsw, model.gsw_max)
    return gsw_out[]
end

# ============================================
# Diagnostics: traits → (fww, s*, s_w) (paper)
# ============================================

"""
    owus_shape_from_climaland(; E0, kx_max, LAI, hc, Td,
                               ks_sat, RAI, dr, Zr,
                               psi_g50, psi_x50, psi_s_sat, b,
                               rho_w=1000.0, g=9.81)

Compute (fww, s_star, s_w) per Bassiouni–Manzoni–Vico (2023).

Units expected:
- E0 [m day^-1]
- kx_max [kg m^-1 MPa^-1 s^-1] (leaf-specific xylem conductivity)
- LAI [m^2 m^-2], hc [m], Td [s day^-1] (usually 86400)
- ks_sat [m day^-1] (soil saturated hydraulic conductivity)
- RAI [m^2 m^-2], dr [m], Zr [m]
- psi_* [MPa] (negative), b [–]
"""
function owus_shape_from_climaland(; E0, kx_max, LAI, hc, Td,
    ks_sat, RAI, dr, Zr, psi_g50, psi_x50, psi_s_sat, b,
    rho_w=1000.0, g=9.81)

    FT = float(promote_type(typeof(E0), typeof(kx_max), typeof(LAI),
                            typeof(hc), typeof(Td)))

    # Upscale conductances per ground area (paper’s Table 1)
    KP_max = FT(kx_max) * FT(LAI) / FT(hc) * (FT(Td) / FT(rho_w))   # [m day^-1 MPa^-1]
    KSR_max = FT(ks_sat) * sqrt(FT(RAI) / (FT(dr) * FT(Zr))) * FT(1e6) / (FT(rho_w) * FT(g))  # [m day^-1 MPa^-1]

    # Π-groups
    ΠR = abs(FT(psi_g50)) / abs(FT(psi_x50))
    ΠF = FT(E0) / (KP_max * abs(FT(psi_g50)))
    ΠT = (KSR_max * abs(FT(psi_s_sat))) / FT(E0)
    ΠS = abs(FT(psi_g50)) / abs(FT(psi_s_sat))

    # fww (closed form)
    halfΠF = ΠF / 2
    rad = (halfΠF + 1)^2 - 2 * ΠF * ΠR
    rad = max(rad, eps(FT))
    fww = 1 - (1 / (2 * ΠR)) * (1 + halfΠF - sqrt(rad))

    # s(β) helper (Eq. 5)
    function s_of_beta(β)
        β = FT(β)  # cast to your promoted float type
        num = (4 * β * ΠS^2 / ΠT) * ((2 * (1 - β) - β * ΠF) / (1 - (1 - β) * ΠR))
        inner = sqrt(max(1 + num, eps(FT))) - 1
        base = (ΠT / (2 * β * ΠS)) * inner
        return (max(base, eps(FT)))^(-FT(1)/FT(b))
    end

    s_star = clamp(s_of_beta(FT(0.95) * fww), FT(0), FT(1))
    s_w    = clamp(s_of_beta(FT(0.05) * fww), FT(0), FT(1))
    return (fww=FT(fww), s_star=s_star, s_w=s_w)
end

# ==============================
# Friendly constructors / glue
# ==============================

"""
    build_owus_from_traits(; kwargs...) -> OWUSStomatalModel

Thin wrapper around `owus_shape_from_climaland` returning OWUSStomatalModel.
"""
function build_owus_from_traits(; kwargs...)
    pars = owus_shape_from_climaland(; kwargs...)
    return OWUSStomatalModel(fww=pars.fww, s_star=pars.s_star, s_w=pars.s_w, gsw_max=Inf)
end

# ------------------------------
# ClimaLand struct auto-plumbing
# ------------------------------
# We avoid hard dependencies on specific ClimaLand types by soft-getting
# common field names with graceful fallbacks. Adjust the symbol lists
# below to match your branch if needed.

# Get first present field among a list of candidate names; return `default` if none.
@inline function _getfirst(x, names::NTuple{N,Symbol}, default=nothing) where {N}
    for nm in names
        if Base.hasproperty(x, nm)
            val = getproperty(x, nm)
            return val
        end
    end
    return default
end

# Unit helpers
@inline _ms_to_mday(x) = x * 86400.0
@inline _mday_to_ms(x) = x / 86400.0

"""
    build_owus_from_ClimaLand(; canopy_params, soil_params, root_params, met_params=nothing, overrides=NamedTuple())

Derive OWUS parameters from typical ClimaLand parameter structs and return
`OWUSStomatalModel`. You can override/force any quantity by passing it
in `overrides` (e.g., `overrides=(E0=2.5, psi_s_sat=-0.003)`).

We expect (but do not require) the following (case-insensitive by symbol):
- canopy_params: LAI, canopy_height/h_c, kx_max/kx_leaf_spec, psi_g50, psi_x50
- soil_params: Ksat/ksat/ks_sat (m s^-1 or m day^-1), psi_sat, b
- root_params: RAI, dr/fine_root_diameter, Zr/rooting_depth
- met_params: E0 (m s^-1 or m day^-1); otherwise we try canopy_params.E0 or compute from PET if provided
"""
function build_owus_from_ClimaLand(; canopy_params, soil_params, root_params, met_params=nothing, overrides=NamedTuple())
    # Pull raw values with generous field-name guesses
    LAI = something(_getfirst(canopy_params, (:LAI, :lai)),  leaf_nothing_error("LAI"))
    hc  = something(_getfirst(canopy_params, (:canopy_height, :h_c, :hc)), leaf_nothing_error("canopy height"))
    kx  = something(_getfirst(canopy_params, (:kx_max, :kx_leaf_spec, :k_x_max, :kxl_spec)), leaf_nothing_error("kx_max"))
    ψg50 = something(_getfirst(canopy_params, (:psi_g50, :ψg50, :psi_g_50, :psi_g50_MPa)), leaf_nothing_error("psi_g50"))
    ψx50 = something(_getfirst(canopy_params, (:psi_x50, :ψx50, :psi_x_50, :psi_x50_MPa)), leaf_nothing_error("psi_x50"))

    RAI = something(_getfirst(root_params, (:RAI, :root_area_index)), 2.0)   # default modestly
    dr  = something(_getfirst(root_params, (:dr, :fine_root_diameter, :root_diameter)), 3e-4)  # m
    Zr  = something(_getfirst(root_params, (:Zr, :rooting_depth, :root_depth)), 1.0)           # m

    # Soil hydraulics
    ks_any = _getfirst(soil_params, (:Ksat, :ksat, :ks_sat, :Ks_sat, :ks_sat_mps, :k_sat))
    ks_sat_mday = if ks_any === nothing
        1.0      # m day^-1 fallback
    else
        k = float(ks_any)
        # heuristic: if very small (<1e-3), assume m s^-1 and convert
        if k < 1e-3
            _ms_to_mday(k)
        else
            k
        end
    end
    ψs_sat = something(_getfirst(soil_params, (:psi_sat, :ψ_sat, :psi_s_sat, :soil_psi_sat)), -0.005)  # MPa
    b      = something(_getfirst(soil_params, (:b, :campbell_b, :soil_b)), 4.38)

    # Atmosphere / demand
    E0_raw = if hasproperty(overrides, :E0)
        overrides.E0
    else
        e0cand = met_params === nothing ? nothing : _getfirst(met_params, (:E0, :PET, :Ep, :E0_mps, :E0_mday))
        e0alt  = _getfirst(canopy_params, (:E0, :Epot, :E_pot))
        e0cand === nothing ? e0alt : e0cand
    end
    if E0_raw === nothing
        # last-resort gentle PET: ~2.5 mm/day
        E0_mday = 2.5e-3    # m day^-1
    else
        e0 = float(E0_raw)
        E0_mday = e0 < 1e-3 ? _ms_to_mday(e0) : e0   # treat small as m s^-1
    end

    # Allow explicit overrides of any inferred quantity
    E0_mday = get(overrides, :E0, E0_mday)
    LAI     = get(overrides, :LAI, LAI)
    hc      = get(overrides, :hc, hc)
    kx      = get(overrides, :kx_max, kx)
    ψg50    = get(overrides, :psi_g50, ψg50)
    ψx50    = get(overrides, :psi_x50, ψx50)
    ψs_sat  = get(overrides, :psi_s_sat, ψs_sat)
    b       = get(overrides, :b, b)
    RAI     = get(overrides, :RAI, RAI)
    dr      = get(overrides, :dr, dr)
    Zr      = get(overrides, :Zr, Zr)

    pars = owus_shape_from_climaland(; E0=E0_mday, kx_max=kx, LAI=LAI, hc=hc, Td=86400.0,
        ks_sat=ks_sat_mday, RAI=RAI, dr=dr, Zr=Zr,
        psi_g50=ψg50, psi_x50=ψx50, psi_s_sat=ψs_sat, b=b)

    return OWUSStomatalModel(fww=pars.fww, s_star=pars.s_star, s_w=pars.s_w, gsw_max=Inf)
end

# Friendly error
leaf_nothing_error(what) = throw(ArgumentError("build_owus_from_ClimaLand: could not infer $what; pass it via overrides=... or ensure it exists in your params."))

end # module

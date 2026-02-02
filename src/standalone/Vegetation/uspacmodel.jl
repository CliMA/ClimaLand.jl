module uSPACStomata

export uSPACStomatalModel,
       stomatal_conductance!,
       uspac_shape_from_climaland,
       build_uspac_from_traits,
       build_uspac_from_ClimaLand,
       uspac_shape_from_Pi,
       pi_groups_from_calibrated_traits,  # exports the helper for model_interface
       extract_soil_params,            # exports soil param extractor
       isohydry_index,
       # Trait distribution extension (Phase 1)
       TraitQuadrature,
       generate_trait_quadrature,
       integrate_over_traits

# Import spatially-varying rooting depth loader
import ClimaCore
import ClimaLand.Canopy: clm_rooting_depth

# ============================================
# Helper functions for extracting fields from p
# ============================================

# --- soft getters ---
@inline _has(x, f::Symbol) = Base.hasproperty(x, f)
@inline _get(x, f::Symbol, dflt=nothing) = _has(x, f) ? getproperty(x,f) : dflt

# Soil saturation s \in [0,1]
@inline function _saturation_field(p, ::Any, ::Type{FT}) where {FT}
    # pick any canopy Field to clone shape/grid from
    like = p.drivers.T

    # 1) Prefer canopy soil_moisture_stress if it exposes θ info
    if Base.hasproperty(p.canopy, :soil_moisture_stress)
        sms = p.canopy.soil_moisture_stress

        if Base.hasproperty(sms, :θ) &&
           Base.hasproperty(sms, :θ_high) &&
           Base.hasproperty(sms, :θ_low)
            θ    = sms.θ
            θ_hi = sms.θ_high
            θ_lo = sms.θ_low
            return @. clamp((θ - θ_lo) / max(θ_hi - θ_lo, eps(FT)), FT(0), FT(1))
        elseif Base.hasproperty(sms, :s)
            s = sms.s
            return @. clamp(FT(s), FT(0), FT(1))
        elseif Base.hasproperty(sms, :βm)
            βm = sms.βm
            # treat βm as a saturation-like 0..1 proxy
            return @. clamp(FT(βm), FT(0), FT(1))
        end
    end

    # 2) Otherwise, try soil state with θ_r, ν
    if Base.hasproperty(p, :soil)
        soil = p.soil
        if Base.hasproperty(soil, :θ) && Base.hasproperty(soil, :parameters)
            θ = soil.θ
            pars = soil.parameters
            if Base.hasproperty(pars, :ν) && Base.hasproperty(pars, :θ_r)
                ν   = pars.ν
                θ_r = pars.θ_r
                return @. clamp((θ - θ_r) / max(ν - θ_r, eps(FT)), FT(0), FT(1))
            end
        end
    end
end

@inline function _E0_field(p, ::Any, ::Type{FT}) where {FT}
    if Base.hasproperty(p.canopy, :energy) && Base.hasproperty(p.canopy.energy, :E0)
        p.canopy.energy.E0
    elseif Base.hasproperty(p.drivers, :E0)
        p.drivers.E0
    else
        # fill a Field with a small constant default
        LAI = p.canopy.biomass.area_index.leaf
        @. zero(LAI) + FT(2.5e-3/86400)
    end
end

@inline function _VPD_field(p, ::Type{FT}) where {FT}
    T = p.drivers.T
    P = p.drivers.P
    q = p.drivers.q
    @inline _svp_pa(Tk) = FT(610.94) * exp((FT(17.625) * (Tk - FT(273.15))) / (Tk - FT(273.15) + FT(243.04)))
    e  = @. (q * P) / (FT(0.622) + (FT(1) - FT(0.622)) * q)
    es = @. _svp_pa(T)
    @. max(es - e, FT(0))
end

# From Bassiouni et al. 2023

# The uSPAC method is fundamentally different from traditional empirical stomatal conductance models 
# in that its parameters have clear biophysical meaning. uSpac is empirical in the sense that β(s) 
# is a fitted piecewise function, but unlike traditional empirical models, 
# its parameters encode measurable plant hydraulic strategies, making them interpretable, 
# predictable, and transferable across species/climates.

# THIS IS NOT OPTIMALITY THEORY AS CODED HERE, but it CAN replicate optimality when appropriate!
# unitless Soil-Plant-Atmosphere continuum (uSPAC; Bassiouni et al. 2023) uses dimensional analysis (Buckingham Π theorem)
# to express stomatal behavior via dimensionless Π-groups representing ratios of:
#   - Plant hydraulic conductances (xylem, guard cells, roots)
#   - Soil-plant-atmosphere pressure gradients
#   - determine whether water use is “supply-limited” or “demand-limited”
# 
# The piecewise linear β(s) form is an empirical stress function that:
#   1. Collapses well across ecosystems when plotted in Π-space
#   2. Can be diagnosed from measurable plant/soil traits
#   3. Can be calibrated from observations (via EKP)

# Flow chart of stomatal conductance calculation:
# Soil Moisture (θ or ψ_soil from ClimaLand)
#     ↓ (converted to)
# Soil Saturation: s = (θ - θ_r) / (θ_sat - θ_r)
#     ↓ (fed into)
# β(s) = piecewise linear stress function
#     ↓ (scales)
# Transpiration: T = β(s) × E₀ × leaf_area_factor
#     ↓ (converted to)
# Molar flux: E_mol = (ρ_water × T) / M_water
#     ↓ (used to calculate)
# Stomatal Conductance: gsw = E_mol × (P_air / VPD)

# Key Points
# Direct scaling: Soil moisture directly scales transpiration via β(s), not through an intermediate optimization step
# No photosynthesis coupling: Unlike Medlyn, this doesn't involve carbon assimilation (An) in the conductance calculation
# Diagnostic approach: Given s, the model immediately knows the stress level and resulting conductance
# Shape matters: The parameters (fww, s_star, s_w) determine:
#   fww: Maximum fraction of potential ET when well-watered
#   s_star: When does water limitation begin?
#   s_w: When does the plant completely shut down?

# Parameters have biophysical meaning: Each Π-group represents a measurable trait ratio.

# Parameters are predictable: Given plant/soil traits, you can compute (fww, s*, s_w) mechanistically.

# Parameters vary systematically:
    # Isohydric species → high ΠR → lower fww, higher s_star
    # Deep-rooted species → high ΠT → lower s_w (can access deeper water)
    # This matches observed trait spectra (Wright et al. 2004, Bartlett et al. 2016).

# Parameters can acclimate!!!!!!: If climate changes → E0 changes → ΠF changes → fww changes predictably.



"""
    uSPACStomatalModel{FT}

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
struct uSPACStomatalModel{FT}
    fww::FT
    s_star::FT
    s_w::FT
    gsw_max::FT
end

uSPACStomatalModel(; fww::Real=0.6, s_star::Real=0.35, s_w::Real=0.10, gsw_max::Real=Inf) = begin
    FT = promote_type(typeof(fww), typeof(s_star), typeof(s_w), typeof(gsw_max), Float64)
    uSPACStomatalModel{FT}(FT(fww), FT(s_star), FT(s_w), FT(gsw_max))
end

# Physical constants (typed helpers)
const _ρw64 = 1000.0             # kg m^-3
const _Mw64 = 0.01801528         # kg mol^-1
const _R64  = 8.314462618        # J mol^-1 K^-1

@inline _consts(::Type{FT}) where {FT} = (FT(_ρw64), FT(_Mw64), FT(_R64))

# β(s) from paper (Eq. 3): 0, linear, plateau at fww
@inline function _beta_piecewise(s::FT, m::uSPACStomatalModel{FT}) where {FT}
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

Given soil saturation `s` and potential evaporation `E0` [m s^-1],
computes stomatal conductance (mol m^-2 s^-1).
"""
function stomatal_conductance!(
    gsw_out::Base.RefValue{FT},
    model::uSPACStomatalModel{FT};
    s::Real, E0::Real, VPD::Real, P_air::Real, T_air::Real, leaf_factor::Real=1
) where {FT}

    sFT, E0FT, VPDFT, P_airFT, T_airFT, leafFT =
        FT(s), FT(E0), FT(VPD), FT(P_air), FT(T_air), FT(leaf_factor)
    ρw, Mw, _ = _consts(FT)

    β = _beta_piecewise(sFT, model)
    T_mps = β * E0FT * leafFT
    E_mol = (ρw * T_mps) / Mw
    gsw   = VPDFT > eps(FT) ? E_mol * (P_airFT / VPDFT) : zero(FT)
    gsw_out[] = min(gsw, model.gsw_max)
    return gsw_out[]
end

# ============================================
# Π-groups to piecewise linear shape parameters
# ============================================

"""
    uspac_shape_from_Pi(; ΠR, ΠF, ΠT, ΠS, b=4.38, β_star_frac=0.95, β_w_frac=0.05)

Compute (fww, s_star, s_w) directly from Π-groups and soil exponent `b`.
"""
function uspac_shape_from_Pi(; ΠR::Real, ΠF::Real, ΠT::Real, ΠS::Real,
                              b::Real=4.38, β_star_frac::Real=0.95, β_w_frac::Real=0.05)

    FT = float(promote_type(typeof(ΠR), typeof(ΠF), typeof(ΠT), typeof(ΠS),
                            typeof(b), typeof(β_star_frac), typeof(β_w_frac)))

    ΠR, ΠF, ΠT, ΠS = FT(ΠR), FT(ΠF), FT(ΠT), FT(ΠS)
    b = FT(b); β_star_frac = FT(β_star_frac); β_w_frac = FT(β_w_frac)

    # fww (closed form)
    halfΠF = ΠF / 2
    rad = (halfΠF + 1)^2 - 2 * ΠF * ΠR
    rad = max(rad, eps(FT))
    fww = one(FT) - (one(FT) / (2 * ΠR)) * (one(FT) + halfΠF - sqrt(rad))

    # s(β) helper (Eq. 5)
    @inline function s_of_beta(β_in)
        β = clamp(FT(β_in), zero(FT), one(FT))
        denom = one(FT) - (one(FT) - β) * ΠR
        denom = ifelse(abs(denom) < sqrt(eps(FT)), sign(denom) * sqrt(eps(FT)), denom)
        termA = (4 * β * ΠS * ΠS) / max(ΠT, eps(FT))
        termB = (2 * (one(FT) - β) - β * ΠF) / denom
        inner = sqrt(max(one(FT) + termA * termB, eps(FT))) - one(FT)
        base  = (ΠT / (2 * β * ΠS)) * inner
        s = (max(base, eps(FT)))^(-one(FT)/b)
        return clamp(s, zero(FT), one(FT))
    end

    s_star = s_of_beta(β_star_frac * fww)
    s_w    = s_of_beta(β_w_frac   * fww)

    return (fww=FT(fww), s_star=FT(s_star), s_w=FT(s_w))
end

# ============================================
# Soil parameter extraction helper
# ============================================

"""
    extract_soil_params(toml_dict, FT)

Extract soil hydraulic parameters from TOML/parameter dict.
Returns NamedTuple (k_sat, ψ_sat, b, ν, θ_r).
"""
function extract_soil_params(toml_dict, FT)
    # Direct access pattern - toml_dict should support string indexing
    # Use try-catch with fallback to defaults
    function safe_get(name, default)
        try
            val = toml_dict[name]
            return FT(val)
        catch
            return FT(default)
        end
    end
    
    # Saturated hydraulic conductivity (m/s)
    k_sat = safe_get("K_sat", 1e-5)
    if k_sat == FT(1e-5)  # Try alternate names if default was used
        k_sat = safe_get("Ksat", safe_get("ksat", safe_get("ks_sat", 1e-5)))
    end
    
    # Matric potential at saturation (MPa, negative)
    ψ_sat = safe_get("ψ_sat", -0.005)
    if ψ_sat == FT(-0.005)
        ψ_sat = safe_get("psi_sat", safe_get("ψsat", -0.005))
    end
    
    # Campbell exponent
    b = safe_get("b", 4.38)
    if b == FT(4.38)
        b = safe_get("campbell_b", safe_get("clapp_hornberger_b", 4.38))
    end
    
    # Porosity
    ν = safe_get("ν", 0.5)
    if ν == FT(0.5)
        ν = safe_get("nu", safe_get("porosity", safe_get("theta_sat", 0.5)))
    end
    
    # Residual water content
    θ_r = safe_get("θ_r", 0.05)
    if θ_r == FT(0.05)
        θ_r = safe_get("theta_r", safe_get("theta_res", 0.05))
    end
    
    return (k_sat = k_sat, ψ_sat = ψ_sat, b = b, ν = ν, θ_r = θ_r)
end

# ============================================
# Build Π groups from calibrated traits and ClimaLand model components
# ============================================

"""
    pi_groups_from_calibrated_traits(model, FT)

Build dimensionless Π-groups (ΠR, ΠF, ΠT, ΠS) from model components.
"""
function pi_groups_from_calibrated_traits(
    p,
    canopy,
    pars,  # uSPACPiParameters
    kx, ψx50, ΠR,
    FT::Type{<:AbstractFloat}
)
    # Extract soil parameters from p.soil at runtime
    if Base.hasproperty(p, :soil) && Base.hasproperty(p.soil, :params)
        soil_params = p.soil.params
        # Try to extract K_sat (saturated hydraulic conductivity)
        k_sat_mps = if Base.hasproperty(soil_params, :K_sat)
            FT(soil_params.K_sat)
        elseif Base.hasproperty(soil_params, :Ksat)
            FT(soil_params.Ksat)
        else
            FT(1e-5)  # fallback m/s
        end
        
        # Try to extract ψ_sat (soil water potential at saturation)
        ψs_sat = if Base.hasproperty(soil_params, :ψ_sat)
            FT(soil_params.ψ_sat) / FT(1e6)  # Convert Pa to MPa
        elseif Base.hasproperty(soil_params, :psi_sat)
            FT(soil_params.psi_sat) / FT(1e6)
        else
            FT(-0.001)  # fallback MPa
        end
    else
        # Fallback if no soil available
        k_sat_mps = FT(1e-5)  # m/s
        ψs_sat = FT(-0.001)   # MPa
    end
    
    k_sat_mday = k_sat_mps * FT(86400)  # Convert m/s to m/day
    
    # Extract root parameters from hydraulics
    hydraulics = canopy.hydraulics
    RAI_val = if Base.hasproperty(hydraulics.parameters, :RAI)
        FT(hydraulics.parameters.RAI)
    else
        FT(4.0)  # default
    end
    
    Zr_val = if Base.hasproperty(hydraulics.parameters, :rooting_depth)
        FT(hydraulics.parameters.rooting_depth)
    else
        FT(1.0)  # default m
    end
    
    # Compute soil-root supply conductance
    K_SR_max = k_sat_mday * RAI_val / Zr_val
    
    # Get canopy hydraulics parameters
    hydraulics = canopy.hydraulics

    # Get canopy structural parameters
    # Extract canopy height from hydraulics (m)
    hc = FT(hydraulics.compartment_surfaces[end])
    
    # Extract actual LAI from auxiliary state (Field)
    LAI = p.canopy.biomass.area_index.leaf
    
    # Extract reference ET from auxiliary state if available, else use default
    E0_mday = _E0_field(p, canopy, FT)
    
    # Sapwood thickness (m) - typical default value
    Td = FT(0.01)
    
    # Physical constants
    rho_w = FT(1000.0)  # kg/m³
    
    # Compute upscaled xylem conductance (m/day)
    # kx is field-like from calibration, broadcast over it
    KP_max = @. kx * LAI / hc * (Td / rho_w)
    
    # Stomatal water potential at 50% closure (MPa)
    # ψx50 and ΠR are field-like from calibration
    ψg50 = @. ΠR * ψx50
    
    # Safeguard denominators to prevent NaN from division by zero
    # This is critical for extreme conditions (hyperarid, very low LAI, etc.)
    KP_max_safe = @. max(abs(KP_max), FT(1e-10))  # Minimum xylem conductance
    ψg50_safe = @. max(abs(ψg50), FT(1e-6))       # Minimum guard cell pressure
    E0_safe = @. max(abs(E0_mday), FT(1e-8))      # Minimum evaporative demand
    ψs_sat_safe = @. max(abs(ψs_sat), FT(1e-6))   # Minimum soil saturation pressure
    
    # ΠF: Plant water flux control (dimensionless)
    ΠF = @. E0_safe / (KP_max_safe * ψg50_safe)
    
    # ΠT: Supply capacity (dimensionless)
    # = (max supply) / (atmospheric demand)
    ΠT = @. (K_SR_max * ψs_sat_safe) / E0_safe
    
    # ΠS: Soil suitability (dimensionless)  
    # = (guard cell failure pressure) / (soil saturation pressure)
    ΠS = @. ψg50_safe / ψs_sat_safe
    
    return (ΠF, ΠT, ΠS)
end


# """
#     build_uspac_from_Pi(; ΠR, ΠF, ΠT, ΠS, b=4.38, β_star_frac=0.95, β_w_frac=0.05, gsw_max=Inf)

# Convenience builder for direct-Π calibration path.
# """
# function build_uspac_from_Pi(; ΠR::Real, ΠF::Real, ΠT::Real, ΠS::Real,
#                               b::Real=4.38, β_star_frac::Real=0.95, β_w_frac::Real=0.05,
#                               gsw_max::Real=Inf)
#     pars = uspac_shape_from_Pi(; ΠR=ΠR, ΠF=ΠF, ΠT=ΠT, ΠS=ΠS,
#                                 b=b, β_star_frac=β_star_frac, β_w_frac=β_w_frac)
#     return (fww=pars.fww, s_star=pars.s_star, s_w=pars.s_w, gsw_max=gsw_max)
# end

# # Builder for state-based construction (needs uSPACConductancePi type from stomatalconductance.jl)
# build_uspac_from_Pi_from_state(; p, canopy, ΓR, ΓF, ΓT, ΓS,
#     b=4.38, β_star_frac=0.95, β_w_frac=0.05, gsw_max=Inf) =
#     uSPACConductancePi{eltype(p.drivers.P)}(
#         uSPACPiParameters{eltype(p.drivers.P)}(;
#             ΓR=ΓR, ΓF=ΓF, ΓT=ΓT, ΓS=ΓS,
#             b=eltype(p.drivers.P)(b),
#             β_star_frac=eltype(p.drivers.P)(β_star_frac),
#             β_w_frac=eltype(p.drivers.P)(β_w_frac),
#             gsw_max=eltype(p.drivers.P)(gsw_max),
#         )
#     )

# ============================================
# Original trait-based interface
# ============================================

"""
    uspac_shape_from_climaland(; E0, kx_max, LAI, hc, Td, ks_sat, RAI, dr, Zr,
                               psi_g50, psi_x50, psi_s_sat, b, rho_w=1000.0, g=9.81)

Compute (fww, s_star, s_w) per Bassiouni–Manzoni–Vico (2023) from raw traits.
"""
function uspac_shape_from_climaland(; E0, kx_max, LAI, hc, Td,
    ks_sat, RAI, dr, Zr, psi_g50, psi_x50, psi_s_sat, b,
    rho_w=1000.0, g=9.81)

    FT = float(promote_type(typeof(E0), typeof(kx_max), typeof(LAI),
                            typeof(hc), typeof(Td)))

    # Upscale conductances
    KP_max = FT(kx_max) * FT(LAI) / FT(hc) * (FT(Td) / FT(rho_w))
    KSR_max = FT(ks_sat) * sqrt(FT(RAI) / (FT(dr) * FT(Zr))) * FT(1e6) / (FT(rho_w) * FT(g))

    # Π-groups
    ΠR = abs(FT(psi_g50)) / abs(FT(psi_x50))
    ΠF = FT(E0) / (KP_max * abs(FT(psi_g50)))
    ΠT = (KSR_max * abs(FT(psi_s_sat))) / FT(E0)
    ΠS = abs(FT(psi_g50)) / abs(FT(psi_s_sat))

    # fww
    halfΠF = ΠF / 2
    rad = (halfΠF + 1)^2 - 2 * ΠF * ΠR
    rad = max(rad, eps(FT))
    fww = 1 - (1 / (2 * ΠR)) * (1 + halfΠF - sqrt(rad))

    # s(β) helper
    function s_of_beta(β)
        β = FT(β)
        num = (4 * β * ΠS^2 / ΠT) * ((2 * (1 - β) - β * ΠF) / (1 - (1 - β) * ΠR))
        inner = sqrt(max(1 + num, eps(FT))) - 1
        base = (ΠT / (2 * β * ΠS)) * inner
        return (max(base, eps(FT)))^(-FT(1)/FT(b))
    end

    s_star = clamp(s_of_beta(FT(0.95) * fww), FT(0), FT(1))
    s_w    = clamp(s_of_beta(FT(0.05) * fww), FT(0), FT(1))
    return (fww=FT(fww), s_star=s_star, s_w=s_w)
end

"""
    build_uspac_from_traits(; kwargs...) -> uSPACStomatalModel

Thin wrapper around `uspac_shape_from_climaland` returning uSPACStomatalModel.
"""
function build_uspac_from_traits(; kwargs...)
    pars = uspac_shape_from_climaland(; kwargs...)
    return uSPACStomatalModel(fww=pars.fww, s_star=pars.s_star, s_w=pars.s_w, gsw_max=Inf)
end

# ============================================
# ClimaLand struct-based builder
# ============================================

@inline function _getfirst(x, names::NTuple{N,Symbol}, default=nothing) where {N}
    for nm in names
        if Base.hasproperty(x, nm)
            val = getproperty(x, nm)
            return val
        end
    end
    return default
end

@inline _ms_to_mday(x) = x * 86400.0
@inline _mday_to_ms(x) = x / 86400.0

"""
    build_uspac_from_ClimaLand(; canopy_params, soil_params, root_params, met_params=nothing, overrides=NamedTuple())

Derive uSPAC parameters from typical ClimaLand parameter structs.
"""
function build_uspac_from_ClimaLand(; canopy_params, soil_params, root_params, met_params=nothing, overrides=NamedTuple())
    @inline function _get_required(x, names::NTuple{N,Symbol}, human::AbstractString) where {N}
        val = _getfirst(x, names, nothing)
        val === nothing && leaf_nothing_error(human)
        return val
    end

    LAI  = _get_required(canopy_params, (:LAI, :lai), "LAI")
    hc   = _get_required(canopy_params, (:canopy_height, :h_c, :hc), "canopy height")
    kx   = _get_required(canopy_params, (:kx_max, :kx_leaf_spec, :k_x_max, :kxl_spec), "kx_max")
    ψg50 = _get_required(canopy_params, (:psi_g50, :ψg50, :psi_g_50, :psi_g50_MPa), "psi_g50")
    ψx50 = _get_required(canopy_params, (:psi_x50, :ψx50, :psi_x_50, :psi_x50_MPa), "psi_x50")

    RAI = something(_getfirst(root_params, (:RAI, :root_area_index)), 2.0)
    dr  = something(_getfirst(root_params, (:dr, :fine_root_diameter, :root_diameter)), 3e-4)
    Zr  = something(_getfirst(root_params, (:Zr, :rooting_depth, :root_depth)), 1.0)

    ks_any = _getfirst(soil_params, (:Ksat, :ksat, :ks_sat, :Ks_sat, :ks_sat_mps, :k_sat))
    ks_sat_mday = if ks_any === nothing
        1.0
    else
        k = float(ks_any)
        k < 1e-3 ? _ms_to_mday(k) : k
    end
    ψs_sat = something(_getfirst(soil_params, (:psi_sat, :ψ_sat, :psi_s_sat, :soil_psi_sat)), -0.005)
    b      = something(_getfirst(soil_params, (:b, :campbell_b, :soil_b)), 4.38)

    E0_raw = if hasproperty(overrides, :E0)
        overrides.E0
    else
        e0cand = met_params === nothing ? nothing : _getfirst(met_params, (:E0, :PET, :Ep, :E0_mps, :E0_mday))
        e0alt  = _getfirst(canopy_params, (:E0, :Epot, :E_pot))
        e0cand === nothing ? e0alt : e0cand
    end
    if E0_raw === nothing
        E0_mday = 2.5e-3
    else
        e0 = float(E0_raw)
        E0_mday = e0 < 1e-3 ? _ms_to_mday(e0) : e0
    end

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

    pars = uspac_shape_from_climaland(; E0=E0_mday, kx_max=kx, LAI=LAI, hc=hc, Td=86400.0,
        ks_sat=ks_sat_mday, RAI=RAI, dr=dr, Zr=Zr,
        psi_g50=ψg50, psi_x50=ψx50, psi_s_sat=ψs_sat, b=b)

    return uSPACStomatalModel(fww=pars.fww, s_star=pars.s_star, s_w=pars.s_w, gsw_max=Inf)
end

leaf_nothing_error(what) = throw(ArgumentError("build_uspac_from_ClimaLand: could not infer $what; pass it via overrides=... or ensure it exists in your params."))

# ============================================
# Isohydry diagnostic
# ============================================

@inline function _ols_slope_intercept_R2(x::AbstractVector{T}, y::AbstractVector{T}) where {T<:Real}
    n = length(x)
    n == length(y) || throw(ArgumentError("x and y must have same length"))
    n ≥ 2 || throw(ArgumentError("need ≥2 points to fit a slope"))

    mx = sum(x) / n
    my = sum(y) / n

    sxx = zero(T); sxy = zero(T); syy = zero(T)
    @inbounds for i in eachindex(x)
        dx = x[i] - mx
        dy = y[i] - my
        sxx += dx*dx
        sxy += dx*dy
        syy += dy*dy
    end

    sxx == 0 && throw(ArgumentError("x has zero variance; slope undefined"))
    m  = sxy / sxx
    a  = my - m*mx
    R2 = syy == 0 ? one(T) : (sxy*sxy)/(sxx*syy)
    return (slope = T(m), intercept = T(a), R2 = T(R2), n = n)
end

@inline function _psi_leaf_linear(E_mol::T, ψs::T; K_lin::T) where {T<:Real}
    return ψs - E_mol / max(K_lin, eps(T))
end

function _psi_leaf_weibull(E_mol::T, ψs::T; Kmax::T, P50::T, a::T, ψ_lo::T=-15.0, ψ_hi::T=0.0) where {T<:Real}
    lo = T(min(ψ_lo, ψs - T(1e-6)))
    hi = T(min(ψ_hi, ψs - T(1e-9)))
    K_of(ψ) = Kmax / (one(T) + exp(-a * (ψ - P50)))
    F(ψ) = K_of(ψ) * (ψs - ψ) - E_mol

    Flo = F(lo); Fhi = F(hi)
    if Flo*Fhi > 0
        return _psi_leaf_linear(E_mol, ψs; K_lin = Kmax)
    end

    max_iter = 60
    for _ in 1:max_iter
        mid = (lo + hi)/2
        Fm  = F(mid)
        if abs(Fm) < eps(T)*max(one(T), abs(E_mol))
            return mid
        end
        if Flo*Fm <= 0
            hi = mid; Fhi = Fm
        else
            lo = mid; Flo = Fm
        end
    end
    return (lo + hi)/2
end

"""
    isohydry_index(model::uSPACStomatalModel{FT}; s, E0, VPD, P_air, T_air, ψ_soil, ...)
        -> (slope, intercept, R2, n)

Compute isohydry slope dΨ_leaf/dΨ_soil.
"""
function isohydry_index(
    model::uSPACStomatalModel{FT};
    s::AbstractVector,
    E0::AbstractVector,
    VPD::AbstractVector,
    P_air::AbstractVector,
    T_air::AbstractVector,
    ψ_soil::AbstractVector,
    leaf_factor::Real = 1,
    hydraulics::Symbol = :linear,
    K_lin::Real = NaN,
    Kmax::Real = NaN, P50::Real = NaN, a::Real = NaN,
) where {FT}

    n = length(s)
    @assert length(E0)==n && length(VPD)==n && length(P_air)==n && length(T_air)==n && length(ψ_soil)==n

    ρw, Mw, _ = _consts(FT)
    Ψl = Vector{FT}(undef, n)
    Ψs = FT.(ψ_soil)

    gref = Ref(FT(0))
    @inbounds for i in 1:n
        s_i    = FT(s[i]);  E0_i = FT(E0[i]);  VPD_i = FT(VPD[i])
        P_i    = FT(P_air[i]);  T_i = FT(T_air[i])
        stomatal_conductance!(gref, model; s=s_i, E0=E0_i, VPD=VPD_i, P_air=P_i, T_air=T_i, leaf_factor=leaf_factor)
        gsw_i  = gref[]
        E_mol  = gsw_i * (VPD_i / P_i)

        if hydraulics === :linear
            isfinite(K_lin) || throw(ArgumentError("K_lin must be provided"))
            Ψl[i] = _psi_leaf_linear(E_mol, Ψs[i]; K_lin=FT(K_lin))
        elseif hydraulics === :weibull
            (isfinite(Kmax) && isfinite(P50) && isfinite(a)) ||
                throw(ArgumentError("Kmax, P50, a must be provided"))
            Ψl[i] = _psi_leaf_weibull(E_mol, Ψs[i]; Kmax=FT(Kmax), P50=FT(P50), a=FT(a))
        else
            throw(ArgumentError("hydraulics must be :linear or :weibull"))
        end
    end

    return _ols_slope_intercept_R2(Ψs, Ψl)
end

# ============================================
# Trait Distribution Extension (Phase 1)
# ============================================

"""
    TraitQuadrature{FT, N}

Quadrature representation of trait distribution for efficient integration.
Stores N sample points in 3D trait space [log(kx), P50, ΠR] with weights.

# Fields
- `points::NTuple{N, NTuple{3, FT}}`: Trait samples (log(kx), P50, ΠR)
- `weights::NTuple{N, FT}`: Quadrature weights (sum to 1)

# Example
```julia
quad = generate_trait_quadrature(
    μ_logkx=0.7, μ_P50=-2.0, μ_ΠR=0.5,
    σ_logkx=0.5, σ_P50=1.5, σ_ΠR=0.15,
    ρ_kx_P50=0.7, n_quad=5
)
# Returns ~5-15 strategically placed points in 3D trait space
```
"""
struct TraitQuadrature{FT <: AbstractFloat, N}
    points::NTuple{N, NTuple{3, FT}}   # (log(kx), P50, ΠR) samples
    weights::NTuple{N, FT}              # Integration weights
end

"""
    generate_trait_quadrature(; μ_logkx, μ_P50, μ_ΠR, σ_logkx, σ_P50, σ_ΠR,
                               ρ_kx_P50=0.0, ρ_kx_ΠR=0.0, ρ_P50_ΠR=0.0,
                               n_quad=5, FT=Float64)

Generate Gaussian-Hermite quadrature points for 3D trait distribution.

# Arguments
- `μ_logkx, μ_P50, μ_ΠR`: Mean traits
- `σ_logkx, σ_P50, σ_ΠR`: Standard deviations  
- `ρ_kx_P50, ρ_kx_ΠR, ρ_P50_ΠR`: Trait correlations (default 0 = independent)
- `n_quad`: Quadrature order (3, 5, or 7 recommended)

# Returns
`TraitQuadrature{FT, N}` with N ≈ n_quad points optimally placed in trait space.

# Physical meaning of correlations
- `ρ_kx_P50 > 0`: Safety-efficiency tradeoff (high kx ↔ less negative P50)
- `ρ_kx_ΠR < 0`: High conductance requires tight regulation (high kx ↔ low ΠR)
- `ρ_P50_ΠR < 0`: Resistant plants can afford to be anisohydric
"""
function generate_trait_quadrature(;
    μ_logkx::Real,
    μ_P50::Real,
    μ_ΠR::Real,
    σ_logkx::Real,
    σ_P50::Real,
    σ_ΠR::Real,
    ρ_kx_P50::Real = 0.0,
    ρ_kx_ΠR::Real = 0.0,
    ρ_P50_ΠR::Real = 0.0,
    n_quad::Int = 5,
    FT::Type = Float64
)
    # Get 1D Gauss-Hermite nodes for standard normal
    x_nodes, w_nodes = _hermite_nodes_weights(n_quad, FT)
    
    # Build covariance matrix from correlations
    ρ_matrix = FT[1.0      ρ_kx_P50  ρ_kx_ΠR;
                  ρ_kx_P50  1.0       ρ_P50_ΠR;
                  ρ_kx_ΠR   ρ_P50_ΠR  1.0]
    
    σ_diag = FT[σ_logkx, σ_P50, σ_ΠR]
    
    # Cholesky decomposition for transform: Σ = L*L^T
    # Build covariance Σ = D*R*D where D=diag(σ), R=correlation matrix
    Σ = zeros(FT, 3, 3)
    for i in 1:3, j in 1:3
        Σ[i,j] = σ_diag[i] * ρ_matrix[i,j] * σ_diag[j]
    end
    
    # Add regularization for numerical stability
    for i in 1:3
        Σ[i,i] += FT(1e-10)
    end
    
    # Simple strategy: Central point + axis-aligned perturbations
    # (Full tensor product would be n^3 points - too many)
    
    μ_vec = FT[μ_logkx, μ_P50, μ_ΠR]
    
    points = NTuple{3, FT}[]
    weights = FT[]
    
    # Central point (mean) - highest weight
    push!(points, (μ_logkx, μ_P50, μ_ΠR))
    push!(weights, w_nodes[div(n_quad, 2) + 1])
    
    # Points along each principal axis (±σ)
    for dim in 1:3
        for i in 1:n_quad
            if i != div(n_quad, 2) + 1  # Skip center
                δ = zeros(FT, 3)
                δ[dim] = sqrt(FT(2)) * σ_diag[dim] * x_nodes[i]
                
                trait_sample = μ_vec + δ
                
                # Clamp to physical bounds
                trait_clamped = _clamp_trait_sample(trait_sample, FT)
                
                push!(points, (trait_clamped[1], trait_clamped[2], trait_clamped[3]))
                push!(weights, w_nodes[i] / FT(3))  # Distribute weight across axes
            end
        end
    end
    
    # Normalize weights to sum to 1
    w_sum = sum(weights)
    weights_normalized = [w / w_sum for w in weights]
    
    N = length(points)
    return TraitQuadrature{FT, N}(
        ntuple(i -> points[i], N),
        ntuple(i -> weights_normalized[i], N)
    )
end

"""
    _hermite_nodes_weights(n::Int, ::Type{FT})

Get 1D Gauss-Hermite quadrature nodes and weights for N(0,1).
Hardcoded values for n=3,5,7 (most efficient orders).
"""
function _hermite_nodes_weights(n::Int, ::Type{FT}) where {FT}
    if n == 3
        x = FT[-1.2247448713915889, 0.0, 1.2247448713915889]
        w = FT[0.29540897515091936, 0.4091820497081613, 0.29540897515091936]
    elseif n == 5
        x = FT[-2.0201828704560856, -0.9585724646138185, 0.0, 
               0.9585724646138185, 2.0201828704560856]
        w = FT[0.019953242059045913, 0.39361932315224116, 0.17245533921312615,
               0.39361932315224116, 0.019953242059045913]
    elseif n == 7
        x = FT[-2.6519613568352334, -1.6735516287674714, -0.8162878828589647, 0.0,
               0.8162878828589647, 1.6735516287674714, 2.6519613568352334]
        w = FT[0.0009717812450995193, 0.05451558281912703, 0.42560725261010205,
               0.03810862634722653, 0.42560725261010205, 0.05451558281912703,
               0.0009717812450995193]
    else
        @warn "Unsupported quadrature order $n, using n=3"
        return _hermite_nodes_weights(3, FT)
    end
    return (x, w)
end

"""
    _clamp_trait_sample(trait::AbstractVector{FT}, ::Type{FT})

Enforce physical bounds on trait samples:
- log(kx): No bounds (any real number valid)
- P50: Must be negative (< 0 MPa)
- ΠR: Must be in [0, 1]
"""
@inline function _clamp_trait_sample(trait::AbstractVector{FT}, ::Type{FT}) where {FT}
    logkx = trait[1]  # No clamping needed
    P50 = min(trait[2], -FT(0.01))  # Must be negative (at least -0.01 MPa)
    ΠR = clamp(trait[3], FT(0), FT(1))  # Bounded [0,1]
    return FT[logkx, P50, ΠR]
end

"""
    integrate_over_traits(f::Function, quad::TraitQuadrature{FT, N})

Compute weighted integral over trait distribution:
```
E[f(traits)] ≈ ∑ᵢ wᵢ × f(trait_sampleᵢ)
```

# Arguments
- `f`: Function taking `(kx, P50, ΠR)` and returning scalar or Field
- `quad`: TraitQuadrature with sample points and weights

# Example
```julia
# Compute ecosystem conductance as weighted average over trait types
quad = generate_trait_quadrature(μ_logkx=0.7, μ_P50=-2.0, μ_ΠR=0.5, ...)

gs_ecosystem = integrate_over_traits(quad) do kx, P50, ΠR
    # Your existing mechanistic uSPAC calculation for this trait type
    ΠF, ΠT, ΠS = compute_pi_groups(kx, P50, ΠR, E0, ...)
    fww, sstar, sw = uspac_shape_from_Pi(ΠR=ΠR, ΠF=ΠF, ΠT=ΠT, ΠS=ΠS)
    β = compute_beta(s, fww, sstar, sw)
    return conductance_from_beta(β, E0, VPD, ...)
end
```
"""
function integrate_over_traits(f::Function, quad::TraitQuadrature{FT, N}) where {FT, N}
    result = zero(FT)
    
    for i in 1:N
        # Extract trait sample: (log(kx), P50, ΠR)
        logkx, P50, ΠR = quad.points[i]
        kx = exp(logkx)  # Transform back from log-space
        
        # Evaluate function for this trait type
        f_val = f(kx, P50, ΠR)
        
        # Weighted sum
        result += quad.weights[i] * f_val
    end
    
    return result
end

end # module

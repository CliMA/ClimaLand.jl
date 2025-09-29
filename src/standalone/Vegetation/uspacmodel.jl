module uSPACStomata

export uSPACStomatalModel,
       stomatal_conductance!,
       uspac_shape_from_climaland,
       build_uspac_from_traits,
       build_uspac_from_ClimaLand,
        uspac_shape_from_Pi, 
        build_uspac_from_Pi,
        build_uspac_from_Pi_from_state

# From Bassiouni et al. 2023

# uSPAC frames the problem in terms of dimensionless Π-groups (ratios of soil, xylem, guard cell, and atmospheric conductances/pressures). These ratios determine whether water use is “supply-limited” or “demand-limited.”
# shows that the optimal strategy in this trait-space is equivalent to a piecewise β(s): flat at high moisture (fww), declining to shutoff at s_w. So the “optimal control” solution reduces to something that looks like an empirical stress function, but whose shape is diagnosed from traits.
# "optimality-based closure" is actually operationalized as a trait-linked stress function instead of a per-timestep optimization. Parameters can be tied directly to measurable hydraulic traits

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

Given soil saturation `s` and potential evaporation `E0` [m s^-1] at canopy,
computes the stomatal conductance (mol m^-2 s^-1) that realizes
T = β(s)*E0, via E_mol = (ρw*T)/Mw and E_mol = g_sw * (VPD/P_air).
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
    T_mps = β * E0FT * leafFT                  # [m s^-1]
    E_mol = (ρw * T_mps) / Mw                  # [mol m^-2 s^-1]
    gsw   = VPDFT > eps(FT) ? E_mol * (P_airFT / VPDFT) : zero(FT)
    gsw_out[] = min(gsw, model.gsw_max)
    return gsw_out[]
end

# ============================================
# Π-space interface: calibrate on ΠR, ΠF, ΠT, ΠS
# ============================================

"""
    uspac_shape_from_Pi(; ΠR, ΠF, ΠT, ΠS, b=4.38, β_star_frac=0.95, β_w_frac=0.05)

Compute (fww, s_star, s_w) directly from Π-groups and soil exponent `b`.

Inputs (dimensionless):
- ΠR = |ψ_g50| / |ψ_x50|
- ΠF = E0 / (K_P,max * |ψ_g50|)
- ΠT = (K_SR,max * |ψ_s_sat|) / E0
- ΠS = |ψ_g50| / |ψ_s_sat|

Other:
- `b`        : Campbell/Clapp–Hornberger soil exponent (default 4.38)
- `β_star_frac`: fraction of fww at which `s_star` is defined (default 0.95)
- `β_w_frac`   : fraction of fww at which `s_w`   is defined (default 0.05)

Returns: NamedTuple `(fww, s_star, s_w)` in promoted float type.
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

"""
    build_uspac_from_Pi(; ΠR, ΠF, ΠT, ΠS, b=4.38, β_star_frac=0.95, β_w_frac=0.05, gsw_max=Inf)
        -> uSPACConductanceParameters-compatible (fww,s*,s_w) pack or your model object

Convenience builder for the direct-Π calibration path.
Use this ONLY when your optimizer proposes Π’s directly.
"""
function build_uspac_from_Pi(; ΠR::Real, ΠF::Real, ΠT::Real, ΠS::Real,
                              b::Real=4.38, β_star_frac::Real=0.95, β_w_frac::Real=0.05,
                              gsw_max::Real=Inf)
    pars = uspac_shape_from_Pi(; ΠR=ΠR, ΠF=ΠF, ΠT=ΠT, ΠS=ΠS,
                                b=b, β_star_frac=β_star_frac, β_w_frac=β_w_frac)
    # If you want a full model object here, adapt to your type; otherwise just return the trio:
    return (fww=pars.fww, s_star=pars.s_star, s_w=pars.s_w, gsw_max=gsw_max)
end

# builder that only sets Γ’s:
build_uspac_from_Pi_from_state(; p, canopy, ΓR, ΓF, ΓT, ΓS,
    b=4.38, β_star_frac=0.95, β_w_frac=0.05, gsw_max=Inf) =
    uSPACConductancePi{eltype(p.drivers.P)}(
        uSPACPiParameters{eltype(p.drivers.P)}(;
            ΓR=ΓR, ΓF=ΓF, ΓT=ΓT, ΓS=ΓS,
            b=eltype(p.drivers.P)(b),
            β_star_frac=eltype(p.drivers.P)(β_star_frac),
            β_w_frac=eltype(p.drivers.P)(β_w_frac),
            gsw_max=eltype(p.drivers.P)(gsw_max),
        )
    )
# update_canopy_conductance!(p, Y, model, canopy) auto-fetches aridity & sand.


# ============================================
# Diagnostics: traits → (fww, s*, s_w)
# ============================================

"""
    uspac_shape_from_climaland(; E0, kx_max, LAI, hc, Td,
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
function uspac_shape_from_climaland(; E0, kx_max, LAI, hc, Td,
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
    build_uspac_from_traits(; kwargs...) -> uSPACStomatalModel

Thin wrapper around `uspac_shape_from_climaland` returning uSPACStomatalModel.
"""
function build_uspac_from_traits(; kwargs...)
    pars = uspac_shape_from_climaland(; kwargs...)
    return uSPACStomatalModel(fww=pars.fww, s_star=pars.s_star, s_w=pars.s_w, gsw_max=Inf)
end

# ------------------------------
# ClimaLand struct
# ------------------------------
# avoid hard dependencies on specific ClimaLand types by soft-getting common field names

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
    build_uspac_from_ClimaLand(; canopy_params, soil_params, root_params, met_params=nothing, overrides=NamedTuple())

Derive uSPAC parameters from typical ClimaLand parameter structs and return
`uSPACStomatalModel`. You can override/force any quantity by passing it
in `overrides` (e.g., `overrides=(E0=2.5, psi_s_sat=-0.003)`).

We expect (but do not require) the following (case-insensitive by symbol):
- canopy_params: LAI, canopy_height/h_c, kx_max/kx_leaf_spec, psi_g50, psi_x50
- soil_params: Ksat/ksat/ks_sat (m s^-1 or m day^-1), psi_sat, b
- root_params: RAI, dr/fine_root_diameter, Zr/rooting_depth
- met_params: E0 (m s^-1 or m day^-1); otherwise we try canopy_params.E0 or compute from PET if provided
"""
function build_uspac_from_ClimaLand(; canopy_params, soil_params, root_params, met_params=nothing, overrides=NamedTuple())
    # Pull raw values with generous field-name guesses
    @inline function _get_required(x, names::NTuple{N,Symbol}, human::AbstractString) where {N}
        val = _getfirst(x, names, nothing)
        val === nothing && leaf_nothing_error(human)
        return val
    end

    LAI  = _get_required(canopy_params, (:LAI, :lai), "LAI") #exists
    hc   = _get_required(canopy_params, (:canopy_height, :h_c, :hc), "canopy height")
    kx   = _get_required(canopy_params, (:kx_max, :kx_leaf_spec, :k_x_max, :kxl_spec), "kx_max")
    ψg50 = _get_required(canopy_params, (:psi_g50, :ψg50, :psi_g_50, :psi_g50_MPa), "psi_g50")
    ψx50 = _get_required(canopy_params, (:psi_x50, :ψx50, :psi_x_50, :psi_x50_MPa), "psi_x50")

    RAI = something(_getfirst(root_params, (:RAI, :root_area_index)), 2.0)   # default modestly # exists
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

    pars = uspac_shape_from_climaland(; E0=E0_mday, kx_max=kx, LAI=LAI, hc=hc, Td=86400.0,
        ks_sat=ks_sat_mday, RAI=RAI, dr=dr, Zr=Zr,
        psi_g50=ψg50, psi_x50=ψx50, psi_s_sat=ψs_sat, b=b)

    return uSPACStomatalModel(fww=pars.fww, s_star=pars.s_star, s_w=pars.s_w, gsw_max=Inf)
end

# Friendly error
leaf_nothing_error(what) = throw(ArgumentError("build_uspac_from_ClimaLand: could not infer $what; pass it via overrides=... or ensure it exists in your params."))

# -------- Isohydry diagnostic (Ψ_leaf vs Ψ_soil slope) -----------------------

export isohydry_index

# Simple OLS fit: y = a + m x
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
    R2 = syy == 0 ? one(T) : (sxy*sxy)/(sxx*syy)  # coefficient of determination
    return (slope = T(m), intercept = T(a), R2 = T(R2), n = n)
end

# Solve E = K(Ψl) * (Ψs - Ψl) for Ψl
# Variant A: linear K (constant)
@inline function _psi_leaf_linear(E_mol::T, ψs::T; K_lin::T) where {T<:Real}
    # Ψl = Ψs - E/K
    return ψs - E_mol / max(K_lin, eps(T))
end

# Variant B: Weibull/sigmoidal K(Ψl) = Kmax / (1 + exp(-a*(Ψl - P50)))
# Solve F(Ψl) = K(Ψl)*(Ψs - Ψl) - E = 0 by bisection on [Ψ_lo, Ψ_hi]
function _psi_leaf_weibull(E_mol::T, ψs::T; Kmax::T, P50::T, a::T, ψ_lo::T=-15.0, ψ_hi::T=0.0) where {T<:Real}
    # Ensure bracket makes sense (Ψ_l expected ≤ Ψ_s; both typically ≤ 0)
    lo = T(min(ψ_lo, ψs - T(1e-6)))
    hi = T(min(ψ_hi, ψs - T(1e-9)))  # slightly below ψs to avoid zero drop singularities
    K_of(ψ) = Kmax / (one(T) + exp(-a * (ψ - P50)))
    F(ψ) = K_of(ψ) * (ψs - ψ) - E_mol

    Flo = F(lo); Fhi = F(hi)
    # If not bracketed, try to widen slightly; if still not, just return linear fallback
    if Flo*Fhi > 0
        return _psi_leaf_linear(E_mol, ψs; K_lin = Kmax)  # conservative fallback
    end

    # Bisection
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
    isohydry_index(
        model::uSPACStomatalModel{FT};
        s::AbstractVector,         # soil saturation [0..1]
        E0::AbstractVector,        # potential evaporation [m s^-1]
        VPD::AbstractVector,       # Pa
        P_air::AbstractVector,     # Pa
        T_air::AbstractVector,     # K (not used in this diagnostic, kept for symmetry)
        ψ_soil::AbstractVector,    # MPa (negative)
        leaf_factor::Real = 1,
        hydraulics::Symbol = :linear,   # :linear or :weibull
        K_lin::Real = NaN,              # needed if hydraulics=:linear  (mol m^-2 s^-1 MPa^-1)
        Kmax::Real = NaN, P50::Real = NaN, a::Real = NaN,  # needed if :weibull
    ) -> (slope, intercept, R2, n)

Compute the isohydry slope dΨ_leaf/dΨ_soil by:
1. Using uSPAC to compute transpiration E (molar) from s and E0, VPD, P_air.
2. Inferring Ψ_leaf from a hydraulic closure (linear or Weibull).
3. Fitting Ψ_leaf = intercept + slope * Ψ_soil (OLS).

Returns: NamedTuple (slope, intercept, R2, n).
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
    @assert length(E0)==n && length(VPD)==n && length(P_air)==n && length(T_air)==n && length(ψ_soil)==n "all inputs must have same length"

    ρw, Mw, _ = _consts(FT)

    Ψl = Vector{FT}(undef, n)
    Ψs = FT.(ψ_soil)

    gref = Ref(FT(0))
    @inbounds for i in 1:n
        # 1) uSPAC beta(s) -> transpiration & E (molar)
        s_i    = FT(s[i]);  E0_i = FT(E0[i]);  VPD_i = FT(VPD[i])
        P_i    = FT(P_air[i]);  T_i = FT(T_air[i])
        stomatal_conductance!(gref, model; s=s_i, E0=E0_i, VPD=VPD_i, P_air=P_i, T_air=T_i, leaf_factor=leaf_factor)
        gsw_i  = gref[]
        E_mol  = gsw_i * (VPD_i / P_i)  # mol m^-2 s^-1

        # 2) Hydraulics: infer Ψ_leaf
        if hydraulics === :linear
            isfinite(K_lin) || throw(ArgumentError("K_lin must be provided for hydraulics=:linear"))
            Ψl[i] = _psi_leaf_linear(E_mol, Ψs[i]; K_lin=FT(K_lin))
        elseif hydraulics === :weibull
            (isfinite(Kmax) && isfinite(P50) && isfinite(a)) ||
                throw(ArgumentError("Kmax, P50, a must be provided for hydraulics=:weibull"))
            Ψl[i] = _psi_leaf_weibull(E_mol, Ψs[i]; Kmax=FT(Kmax), P50=FT(P50), a=FT(a))
        else
            throw(ArgumentError("hydraulics must be :linear or :weibull"))
        end
    end

    # 3) Fit Ψ_leaf = intercept + slope * Ψ_soil
    return _ols_slope_intercept_R2(Ψs, Ψl)
end


end # module

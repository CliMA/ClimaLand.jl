export OptimalLAIParameters,
    compute_L_max, compute_m, lambertw0, compute_steady_state_LAI

"""
    OptimalLAIParameters{FT<:AbstractFloat}

The required parameters for the optimal LAI model based on Zhou et al. (2025).

Water limitation is handled through the f0*P/A0 term following Zhou et al. (2025) Equation 11,
where P is annual precipitation and A0 is annual potential GPP.

# References
Zhou et al. (2025) "A General Model for the Seasonal to Decadal Dynamics of Leaf Area"
Global Change Biology. https://onlinelibrary.wiley.com/doi/pdf/10.1111/gcb.70125

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct OptimalLAIParameters{FT <: AbstractFloat}
    """Light extinction coefficient (dimensionless), typically 0.5"""
    k::FT
    """Unit cost of constructing and maintaining leaves (mol m^-2 yr^-1), globally fitted as 12.227 mol m^-2 yr^-1"""
    z::FT
    """Dimensionless parameter representing departure from square-wave LAI dynamics, globally fitted as 0.771"""
    sigma::FT
    """Smoothing factor for exponential moving average (dimensionless, 0-1). Set to 0.067 for ~15 days of memory"""
    alpha::FT
    """Peak fraction of annual precipitation available for transpiration (dimensionless,
    0-1), reached at the energy-water limitation transition. The fraction actually used,
    `f0 = f0_max * exp(-0.604 * ln^2(AI/1.9))` with `AI` the aridity index, falls off
    toward both extremes; Zhou et al. (2025) fit `f0_max = 0.65`. See `f0_from_aridity`."""
    f0_max::FT
    """Long-term memory timescale (s) of the A0 and precipitation running-mean annual
    totals that set LAI_max and the steady-state LAI. Default 2 years; a longer value
    filters the seasonal cycle more strongly, avoiding aliasing of the annual cycle."""
    tau_long_term::FT
    """Whether the C3 fraction is computed from the C3/C4 competition on the
    running-mean per-pathway potential GPP (1.0, the default) or held fixed at the
    photosynthesis model's static value (0.0). See `c3_fraction_from_competition`."""
    online_c3c4::FT
    """Steepness of the logistic mapping the proportional C4 GPP advantage to the
    expected C4 fraction (dimensionless), fitted by Lavergne et al. (2022)."""
    c3c4_k::FT
    """Midpoint of that logistic (dimensionless): the C4 GPP advantage at which C4
    and C3 are equally expected."""
    c3c4_q::FT
    """Coefficient `a` of the C3 tree-cover relation `tc(g) = a·g^b + c`, with `g`
    the annual C3 GPP (kg C m^-2 yr^-1)."""
    tc_a::FT
    """Exponent `b` of the tree-cover relation (dimensionless)."""
    tc_b::FT
    """Offset `c` of the tree-cover relation; negative, so tree cover vanishes below
    a threshold GPP."""
    tc_c::FT
    """Reference annual C3 GPP (kg C m^-2 yr^-1) normalizing the tree-cover relation:
    `tc(g)/tc(tc_gpp_ref)`, clamped to [0, 1], is the C3 tree proportion."""
    tc_gpp_ref::FT
end

Base.eltype(::OptimalLAIParameters{FT}) where {FT} = FT

# make these custom structs broadcastable as tuples
Base.broadcastable(x::OptimalLAIParameters) = tuple(x)

"""
    OptimalLAIParameters{FT}(toml_dict::CP.ParamDict) where {FT}

Creates an `OptimalLAIParameters` object from a TOML parameter dictionary.
"""
function OptimalLAIParameters{FT}(toml_dict::CP.ParamDict) where {FT}
    return OptimalLAIParameters{FT}(
        k = FT(toml_dict["optimal_lai_k"]),
        z = FT(toml_dict["optimal_lai_z"]),
        sigma = FT(toml_dict["optimal_lai_sigma"]),
        alpha = FT(toml_dict["optimal_lai_alpha"]),
        f0_max = FT(toml_dict["optimal_lai_f0_max"]),
        tau_long_term = FT(toml_dict["optimal_lai_tau_long_term"]),
        online_c3c4 = FT(toml_dict["optimal_lai_online_c3c4"]),
        c3c4_k = FT(toml_dict["optimal_lai_c3c4_k"]),
        c3c4_q = FT(toml_dict["optimal_lai_c3c4_q"]),
        tc_a = FT(toml_dict["optimal_lai_tc_a"]),
        tc_b = FT(toml_dict["optimal_lai_tc_b"]),
        tc_c = FT(toml_dict["optimal_lai_tc_c"]),
        tc_gpp_ref = FT(toml_dict["optimal_lai_tc_gpp_ref"]),
    )
end

"""
    compute_L_max(Ao_annual, k, z, precip_annual, f0, ca_pa, chi, vpd_gs)

Compute seasonal maximum leaf area index (LAI_max) based on annual potential GPP
and water availability, following Zhou et al. (2025) Equation 11.

LAI_max is determined by the minimum of energy-limited and water-limited fAPAR:
- Energy-limited: fAPAR_energy = 1 - z/(k*A0)
- Water-limited: fAPAR_water = f0*P/A0 * (ca(1-chi))/(1.6*D)

# Arguments
- `Ao_annual::FT`: Annual total potential GPP (mol CO2 m^-2 yr^-1).
- `k::FT`: Light extinction coefficient (dimensionless), typically 0.5
- `z::FT`: Unit cost of constructing and maintaining leaves (mol m^-2 yr^-1), 12.227
- `precip_annual::FT`: Mean annual precipitation (mol H2O m^-2 yr^-1)
- `f0::FT`: Fraction of precipitation available for transpiration (dimensionless), 0.65
- `ca_pa::FT`: Ambient CO2 partial pressure (Pa), typically ~40 Pa at 400 ppm
- `chi::FT`: Optimal ratio of intercellular to ambient CO2 (dimensionless), typically 0.7-0.8
- `vpd_gs::FT`: Mean vapor pressure deficit during growing season (Pa)

# Returns
- `LAI_max::FT`: Seasonal maximum leaf area index (m^2 m^-2)

# Notes
Following Zhou et al. (2025) Equation 11:
```
fAPAR_max = min{1 - z/(k*A0), f0*P/A0 * (ca(1-chi))/(1.6*D)}
```
The first term is energy-limited (carbon gain vs leaf cost trade-off).
The second term is water-limited (precipitation constrains transpiration, scaled by
intrinsic water use efficiency iWUE = ca(1-chi)/(1.6*D)).

The iWUE factor converts water flux to carbon flux:
- ca(1-chi): CO2 drawdown from ambient to intercellular (Pa)
- 1.6*D: VPD adjusted for CO2/H2O diffusivity ratio (Pa)

# References
Zhou et al. (2025) Global Change Biology, Equation 11
"""
function compute_L_max(
    Ao_annual::FT,      # mol CO2 m^-2 yr^-1
    k::FT,              # dimensionless
    z::FT,              # mol m^-2 yr^-1
    precip_annual::FT,  # mol H2O m^-2 yr^-1
    f0::FT,             # dimensionless
    ca_pa::FT,          # Pa
    chi::FT,            # dimensionless
    vpd_gs::FT,         # Pa
) where {FT}
    # Handle edge case: very small or zero Ao_annual (e.g., polar regions)
    # When Ao_annual ~ 0, z / (k * Ao_annual) -> Inf, causing numerical issues.
    # Use ifelse for GPU compatibility.
    Ao_annual_safe = max(Ao_annual, eps(FT))

    # Energy-limited fAPAR (Equation 11, first term)
    # Plants optimize leaf area to maximize carbon gain minus construction cost
    fAPAR_energy = FT(1) - z / (k * Ao_annual_safe)

    # Water-limited fAPAR (Equation 11, second term)
    # fAPAR_water = f0 * P / A0 * (ca(1-chi)) / (1.6*D)
    # The iWUE factor (ca(1-chi))/(1.6*D) converts water flux to carbon flux
    # Guard against zero VPD
    vpd_safe = max(vpd_gs, eps(FT))
    iWUE_factor = (ca_pa * (FT(1) - chi)) / (FT(1.6) * vpd_safe)
    fAPAR_water = f0 * precip_annual / Ao_annual_safe * iWUE_factor

    # fAPAR_max is the minimum of energy and water constraints (Equation 11)
    fAPAR_max = min(fAPAR_energy, fAPAR_water)

    # Ensure fAPAR is in valid range [0, 1]
    fAPAR_max = max(FT(0), min(FT(1), fAPAR_max))

    # Convert fAPAR to LAI using Beer's law (Equation 12)
    # fAPAR = 1 - exp(-k * LAI)  ->  LAI = -(1/k) * ln(1 - fAPAR)
    # Guard against fAPAR_max = 1 which would give -log(0) = Inf
    fAPAR_max_safe = min(fAPAR_max, FT(1) - eps(FT))
    LAI_max = -(FT(1) / k) * log(FT(1) - fAPAR_max_safe)

    return LAI_max
end

fAPAR_max_fun(k::FT, LAI_max::FT) where {FT} = FT(1) - exp(-k * LAI_max)

"""
    compute_m(GSL, LAI_max, Ao_annual, sigma, k)

Compute the parameter m, which represents the ratio of steady-state LAI to steady-state GPP.

This implements Equation 20 from Zhou et al. (2025). The parameter m quantifies the
relationship between LAI and GPP dynamics, representing the extent to which seasonal LAI
dynamics depart from a "square wave" (where maximum LAI would be maintained throughout
the growing season).

# Arguments
- `GSL::FT`: Growing season length (days). Defined as the length of continuous period
  above 0C longer than 5 days.
- `LAI_max::FT`: Seasonal maximum leaf area index (m^2 m^-2, dimensionless)
- `Ao_annual::FT`: Annual total potential GPP (mol m^-2 yr^-1). This is the integral of daily
  A0 over the year.
- `sigma::FT`: Dimensionless parameter representing departure from square-wave LAI dynamics.
  Globally fitted as sigma = 0.771
- `k::FT`: Light extinction coefficient (dimensionless)

# Returns
- `m::FT`: Parameter relating steady-state LAI to steady-state GPP (dimensionless, units
  work out as: days * m^2 m^-2 / (mol m^-2 yr^-1 * dimensionless) with implicit conversion)

# References
Zhou et al. (2025) Global Change Biology, Equation 20
"""
function compute_m(
    GSL::FT,        # days
    LAI_max::FT,    # m^2 m^-2 (dimensionless)
    Ao_annual::FT,  # mol m^-2 yr^-1
    sigma::FT,      # dimensionless
    k::FT,
) where {FT}
    # Equation 20: m = (sigma * GSL * LAI_max) / (A0_sum * fAPAR_max)
    fAPAR_max = fAPAR_max_fun(k, LAI_max)

    # Guard against division by zero: when LAI_max ~ 0, fAPAR_max ~ 0,
    # but the numerator (sigma * GSL * LAI_max) is also ~ 0, so m ~ 0 naturally.
    fAPAR_max_safe = max(fAPAR_max, eps(FT))
    Ao_annual_safe = max(Ao_annual, eps(FT))
    m = (sigma * GSL * LAI_max) / (Ao_annual_safe * fAPAR_max_safe)
    return m
end

const MINARG = -inv(Base.MathConstants.e)

"""
    _lambertw0_initial_guess(x::T) where {T<:AbstractFloat}

Provide a robust initial guess for the Lambert W0 function for use in iterative solvers.

# Arguments
- `x::T`: Input value, should be >= -1/e

# Returns
- Initial guess for W0(x)

# Algorithm
- For x > 1: uses log(x) - log(log(x)) approximation
- For x < -0.32 (near -1/e): uses series expansion for accurate convergence near branch point
- For -0.32 <= x <= 1: uses max(x, -0.3) as a simple starting point
"""
@inline function _lambertw0_initial_guess(x::T) where {T <: AbstractFloat}
    if x > one(T)
        return log(x) - log(max(log(x), T(1e-6)))
    elseif x < T(-0.32)
        # Near the branch point -1/e, use series expansion
        # This handles the singular behavior at x = -1/e where W(x) = -1
        p = sqrt(T(2) * (T(ℯ) * x + one(T)))
        return -one(T) + p - p^2 / T(3) + p^3 * T(11) / T(72)
    else
        return max(x, T(-0.3))
    end
end

"""
    lambertw0(x::T; maxiter::Int = 8) where {T<:AbstractFloat}

Compute the principal branch (W0) of the Lambert W function for x in [-1/e, Inf).

This is a GPU-device-friendly implementation using a fixed number of Halley iterations.
The Lambert W function satisfies W(x)*exp(W(x)) = x.

# Arguments
- `x::T`: Input value, must be >= -1/e ~ -0.36788
- `maxiter::Int`: Maximum number of Halley iterations (default: 8; Halley's method has cubic convergence, so 8 is generous)

# Returns
- `W::T`: Lambert W0(x), the principal branch value, or NaN for invalid inputs

# Algorithm
Uses Halley's method with a fixed number of iterations for GPU compatibility:
- No dynamic memory allocation
- No conditional breaks (runs all iterations)
- Broadcastable for use with CuArrays: lambertw0.(cuarray)

# Device Compatibility
This implementation is designed to work on both CPU and GPU:
- All operations are scalar and supported on CUDA.jl
- No array allocations or dynamic loops
- Type-generic over AbstractFloat (Float32, Float64)

# References
Corless et al. (1996) "On the Lambert W function"
"""
@inline function lambertw0(x::T; maxiter::Int = 8) where {T <: AbstractFloat}
    # In our usage, arg = -k*mu*exp(-k*mu) with k > 0, mu >= 0,
    # so x is always in [-1/e, 0]. This check is a safety net.
    if !(isfinite(x)) || x < T(MINARG)
        return T(NaN)
    end
    w = _lambertw0_initial_guess(x)
    for i in 1:maxiter
        ew = exp(w)
        f = w * ew - x
        # Halley denominator
        # Special case: when w ~ -1, both numerator and denominator approach 0
        # This happens at the branch point x = -1/e, where W(-1/e) = -1
        w_plus_1 = w + one(T)
        if abs(w_plus_1) < eps(T)
            # Already at or very near the solution w = -1, no update needed
            Δ = zero(T)
        else
            two_w_plus_2 = T(2) * w_plus_1
            if abs(two_w_plus_2) < eps(T)
                # Near w = -1, use Newton's method instead of Halley
                Δ = f / (ew * w_plus_1)
            else
                denom = ew * w_plus_1 - (w + T(2)) * f / two_w_plus_2
                if abs(denom) < eps(T)
                    Δ = f / (ew * w_plus_1)
                else
                    Δ = f / denom
                end
            end
        end
        w -= Δ
    end
    return w
end

"""
    compute_steady_state_LAI(Ao_daily, m, k, LAI_max)

Compute steady-state LAI from daily potential GPP using the Lambert W function solution.

This implements Equations 13-15 from Zhou et al. (2025). The steady-state LAI (L_s) is
the LAI that would be in equilibrium with GPP if weather conditions were held constant.
Given daily meteorological conditions, this is computed on a daily basis.

# Arguments
- `Ao_daily::FT`: Daily potential GPP (mol m^-2 day^-1). This is the GPP that would be
  achieved if fAPAR = 1, calculated from LUE * PPFD.
- `m::FT`: Parameter relating steady-state LAI to steady-state GPP (dimensionless), from
  `compute_m()`
- `k::FT`: Light extinction coefficient (dimensionless), typically 0.5
- `LAI_max::FT`: Seasonal maximum LAI constraint (m^2 m^-2, dimensionless)

# Returns
- `L_steady::FT`: Steady-state leaf area index (m^2 m^-2, dimensionless). Always >= 0.

# Notes
The solution uses the Lambert W0 function: L_s = min{mu + (1/k)W0[-k*mu*exp(-k*mu)], LAI_max}
where mu = m * A0. The result is constrained to be non-negative and below LAI_max.

# References
Zhou et al. (2025) Global Change Biology, Equations 13-15
"""
function compute_steady_state_LAI(
    Ao_daily::FT,  # mol m^-2 day^-1
    m::FT,         # dimensionless
    k::FT,         # dimensionless
    LAI_max::FT,   # m^2 m^-2 (dimensionless)
) where {FT}
    # mu = m * A0 (Equation 15)
    mu = m * Ao_daily

    # Compute argument for Lambert W function
    arg = -k * mu * exp(-k * mu)

    # Check if argument is in valid range for W0 branch: [-1/e, 0]
    # If outside this range, use boundary solution
    if arg < -FT(1) / FT(exp(1))
        # Beyond valid range; use maximum possible LAI
        L_s = LAI_max
    else
        # Equation 15: L_s = mu + (1/k) * W0[-k mu exp(-k mu)]
        # Using our custom lambertw0 function (W0 is the principal branch)
        w_val = lambertw0(arg)
        L_s = mu + (FT(1) / k) * w_val

        # Take minimum with LAI_max (Equation 15)
        L_s = min(L_s, LAI_max)
    end

    # Ensure non-negative (should be guaranteed mathematically, but enforce for numerical stability)
    L_s = max(zero(FT), L_s)

    return L_s
end

"""
    compute_L_steady_target(A0_daily, k, A0_annual, z, GSL, sigma, precip_annual, f0, ca_pa, chi, vpd_gs)

Compute the steady-state LAI target `L_steady` (Zhou et al. 2025 Eqs. 11-15) from
the daily and annual potential GPP and the water-limitation inputs, without the
Eq. 16 acclimation lag. The prognostic `LAI` relaxes toward this target in the
biomass tendency, so the acclimation `alpha` is applied there rather than here.
"""
function compute_L_steady_target(
    A0_daily::FT,
    k::FT,
    A0_annual::FT,
    z::FT,
    GSL::FT,
    sigma::FT,
    precip_annual::FT,
    f0::FT,
    ca_pa::FT,
    chi::FT,
    vpd_gs::FT,
) where {FT}
    LAI_max =
        compute_L_max(A0_annual, k, z, precip_annual, f0, ca_pa, chi, vpd_gs)
    m = compute_m(GSL, LAI_max, A0_annual, sigma, k)
    return compute_steady_state_LAI(A0_daily, m, k, LAI_max)
end

# Shape of Zhou et al. (2025)'s f0(AI) curve: the aridity index at which f0 peaks,
# and the width of the falloff. Structural to the published relation, not tunable.
const AI_PEAK = 1.9
const AI_WIDTH = 0.604

"""
    f0_from_aridity(PET_annual::FT, precip_annual::FT, f0_max::FT) where {FT}

Climate-responsive fraction of precipitation available for transpiration
(Zhou et al. 2025): `f0 = f0_max·exp(−0.604·ln²(AI/1.9))` with aridity index
`AI = PET_annual/precip_annual`. Peaks at `f0_max` at the energy–water transition
(AI = 1.9) and declines toward both the arid and humid extremes.
"""
function f0_from_aridity(
    PET_annual::FT,
    precip_annual::FT,
    f0_max::FT,
) where {FT}
    AI = max(PET_annual, eps(FT)) / max(precip_annual, eps(FT))
    return f0_max * exp(-FT(AI_WIDTH) * log(max(AI, eps(FT)) / FT(AI_PEAK))^2)
end

"""
    aridity_from_f0(f0::FT, f0_max::FT) where {FT}

Inverse of [`f0_from_aridity`](@ref), used to seed `PET_annual` so the online `f0`
starts at the value of the map it replaces. `f0` is symmetric in `ln(AI/1.9)`, so
this returns the arid branch (`AI ≥ 1.9`), which is where the map was fitted. An
`f0` at or above `f0_max` has no arid-branch preimage and returns the peak `1.9`.
"""
function aridity_from_f0(f0::FT, f0_max::FT) where {FT}
    return FT(AI_PEAK) *
           exp(sqrt(max(log(f0_max / max(f0, eps(FT))), FT(0)) / FT(AI_WIDTH)))
end

"""
    c3_fraction_from_competition(A0c3_annual, A0c4_annual, Mc, fapar, parameters)

Online C3 fraction from the pyrealm-style C3/C4 competition (GPP-derived, no PFT),
using the running-mean per-pathway potential GPP `A0c3_annual`/`A0c4_annual`
(mol CO2 m^-2 yr^-1) and the molar mass `Mc` (kg mol^-1). Steps (Lavergne/pyrealm):
the proportional C4 GPP advantage `A4 = (A0c4 − A0c3)/A0c3` is passed through a
logistic to an expected C4 fraction, which is then penalised by the proportion of
C3-tree canopy (estimated from `A0c3`), so C4 is suppressed where C3 trees would
shade it. Returns `fractional_c3 = 1 − frac_c4`.
"""
function c3_fraction_from_competition(
    A0c3_annual::FT,
    A0c4_annual::FT,
    Mc::FT,
    fapar::FT,
    parameters::OptimalLAIParameters{FT},
) where {FT}
    (; c3c4_k, c3c4_q, tc_a, tc_b, tc_c, tc_gpp_ref) = parameters
    a0c3 = max(A0c3_annual, eps(FT))
    # C4 GPP advantage is a ratio, so it is invariant to the potential→actual GPP
    # scaling (fapar cancels); use the potential per-pathway means directly.
    adv = (A0c4_annual - a0c3) / a0c3
    frac_c4 = 1 / (1 + exp(-c3c4_k * (adv / FT(ℯ) - c3c4_q)))
    # C3-tree proportion from annual C3 GPP in kg m^-2 yr^-1. The tc() relation
    # (pyrealm/Lavergne) is fit to REALIZED GPP, so scale the potential a0c3 by the
    # realized fAPAR — using potential GPP here saturates prop_trees and wrongly
    # suppresses C4 in productive grasslands.
    gppc3 = a0c3 * Mc * fapar
    tc(g) = tc_a * g^tc_b + tc_c
    prop_trees = clamp(tc(gppc3) / tc(tc_gpp_ref), FT(0), FT(1))
    frac_c4 *= (1 - prop_trees)
    return 1 - frac_c4
end

"""
    potential_evaporation(SW_d, LW_d, T_air, ϵ_sfc, σ, λv, M_w)

Net-radiation potential evaporation (mol H2O m^-2 s^-1) over a reference vegetated
surface, the numerator of the aridity index `AI = PET_annual/precip_annual` that sets
the climate-responsive `f0`. Uses the net-radiation option of Zhou et al.'s aridity
input, `Rn = (1 - α_ref) SW_d + ϵ_sfc (LW_d - σ T^4)` (W m^-2), clipped at zero so
night-time negative `Rn` does not draw down the running total.

`α_ref` is the FAO-56 reference-crop albedo, which *defines* the reference surface
the `f0(AI)` relation was fitted against; it is deliberately not the simulated
surface albedo, which would decouple `AI` from the curve that consumes it.
"""
function potential_evaporation(
    SW_d::FT,
    LW_d::FT,
    T_air::FT,
    ϵ_sfc::FT,
    σ::FT,
    λv::FT,
    M_w::FT,
) where {FT}
    α_ref = FT(0.23)
    Rn = (1 - α_ref) * SW_d + ϵ_sfc * (LW_d - σ * T_air^4)
    return max(Rn, zero(FT)) / (λv * M_w)
end

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
    for _ in 1:maxiter
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

# Aridity index at which Zhou et al. (2025)'s f0(AI) peaks, and the width of its
# falloff; part of the published relation, not tunable.
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
    return f0_max * exp(-FT(AI_WIDTH) * log(AI / FT(AI_PEAK))^2)
end

"""
    aridity_from_f0(f0::FT, f0_max::FT) where {FT}

Inverse of [`f0_from_aridity`](@ref), used to seed `PET_annual` so the online `f0`
starts at the value of the map it replaces. `f0` is symmetric in `ln(AI/1.9)`, so
this returns the arid branch (`AI ≥ 1.9`), where the map was fitted. An `f0` at or
above `f0_max` has no preimage and returns the peak `1.9`.
"""
function aridity_from_f0(f0::FT, f0_max::FT) where {FT}
    f0 = clamp(f0, eps(FT), f0_max)
    return FT(AI_PEAK) * exp(sqrt(log(f0_max / f0) / FT(AI_WIDTH)))
end

"""
    canopy_composition_from_competition(A0c3_annual, A0c4_annual, GPPc3_annual, Mc, parameters)

Partition of the canopy into C3 trees, C3 grasses and C4 grasses from the C3/C4
competition of Lavergne et al. (2022), as implemented in pyrealm, on the trailing
per-pathway potential GPP `A0c3_annual`/`A0c4_annual` and the trailing realized C3
GPP `GPPc3_annual` (the potential scaled by fAPAR; all mol CO2 m^-2 yr^-1); `Mc` is
the molar mass of carbon (kg mol^-1). Returns a `NamedTuple`
`(; tree, c3_grass, c4_grass)` summing to one.

The fractions are shares of productivity, not of ground area: the proportional C4 GPP
advantage `(A0c4 − A0c3)/A0c3` goes through a logistic to an expected C4 share of
the open canopy, and the tree share is the C3 tree cover estimated from the annual
realized C3 GPP, normalized by the cover at canopy closure (`tc_gpp_ref`). C4 grasses
are shaded out under trees, so the C4 and C3 grass shares are the open-canopy split
scaled by `1 − tree`.
"""
function canopy_composition_from_competition(
    A0c3_annual::FT,
    A0c4_annual::FT,
    GPPc3_annual::FT,
    Mc::FT,
    parameters::OptimalLAIParameters{FT},
) where {FT}
    (; c3c4_k, c3c4_q) = parameters
    a0c3 = max(A0c3_annual, eps(FT))
    adv = (A0c4_annual - a0c3) / a0c3
    # pyrealm scales the advantage by exp(1/(1+TC)) with TC the observed tree
    # cover; with no such input, TC = 0 leaves the divisor ℯ.
    open_c4 = 1 / (1 + exp(-c3c4_k * (adv / FT(ℯ) - c3c4_q)))
    tree = tree_share_from_gpp(GPPc3_annual, Mc, parameters)
    c4_grass = open_c4 * (1 - tree)
    c3_grass = (1 - open_c4) * (1 - tree)
    return (; tree, c3_grass, c4_grass)
end

"""
    c3_fraction_from_competition(A0c3_annual, A0c4_annual, GPPc3_annual, Mc, parameters)

C3 fraction of the canopy, `1 − c4_grass` of
`canopy_composition_from_competition` (trees are all C3).
"""
function c3_fraction_from_competition(
    A0c3_annual::FT,
    A0c4_annual::FT,
    GPPc3_annual::FT,
    Mc::FT,
    parameters::OptimalLAIParameters{FT},
) where {FT}
    composition = canopy_composition_from_competition(
        A0c3_annual,
        A0c4_annual,
        GPPc3_annual,
        Mc,
        parameters,
    )
    return 1 - composition.c4_grass
end

"""
    tree_share_from_gpp(GPPc3_annual, Mc, parameters)

C3 tree share of the canopy in `canopy_composition_from_competition`: the Lavergne
et al. (2022) tree cover `tc(g) = a·g^b + c` at the annual realized C3 GPP
`GPPc3_annual` (mol CO2 m^-2 yr^-1, converted with the molar mass of carbon `Mc`),
relative to the cover at canopy closure `tc_gpp_ref`, clamped to [0, 1].
"""
function tree_share_from_gpp(
    GPPc3_annual::FT,
    Mc::FT,
    parameters::OptimalLAIParameters{FT},
) where {FT}
    (; tc_a, tc_b, tc_c, tc_gpp_ref) = parameters
    # The tree-cover relation is fitted to annual realized GPP in kg C m^-2 yr^-1.
    gppc3 = max(GPPc3_annual, FT(0)) * Mc
    tc(g) = tc_a * g^tc_b + tc_c
    return clamp(tc(gppc3) / tc(tc_gpp_ref), FT(0), FT(1))
end

"""
    c4_advantage_for_c3_fraction(fractional_c3, tree, parameters)

Inverse of the C3/C4 competition: the proportional C4 GPP advantage
`(A0c4 − A0c3)/A0c3` for which `canopy_composition_from_competition` returns the C3
fraction `fractional_c3`, given the tree share `tree`. Used to seed `A0c4_annual` so
the competition starts at a prescribed C3 map.

The advantage is bounded below by −1 (`A0c4 ≥ 0`), so a pure-C3 cell keeps a small
C4 grass share, and a C4 grass share above the open canopy `1 − tree` has no preimage
and saturates.
"""
function c4_advantage_for_c3_fraction(
    fractional_c3::FT,
    tree::FT,
    parameters::OptimalLAIParameters{FT},
) where {FT}
    (; c3c4_k, c3c4_q) = parameters
    δ = sqrt(eps(FT))
    open_c4 = clamp((1 - fractional_c3) / max(1 - tree, δ), δ, 1 - δ)
    adv = FT(ℯ) * (c3c4_q + log(open_c4 / (1 - open_c4)) / c3c4_k)
    return max(adv, -one(FT))
end

# FAO-56 (Allen et al., 1998) reference crop, which defines the PET that f0(AI) was
# fitted with: albedo, bulk surface resistance (s m^-1), the aerodynamic resistance
# numerator (r_a = FAO56_RA_WIND / u_2 in s m^-1), and the log-profile coefficients of
# their Eq. 47 mapping wind at height h to 2 m, u_2 = u·a / ln(b·h − c).
const FAO56_ALBEDO = 0.23
const FAO56_SURFACE_RESISTANCE = 70
const FAO56_RA_WIND = 208
const FAO56_WIND_A = 4.87
const FAO56_WIND_B = 67.8
const FAO56_WIND_C = 5.42

"""
    potential_evaporation(
        SW_d, LW_d, T_air, P_air, q_air, u_air, h_atmos, ϵ_sfc, σ, M_w, thermo_params,
    )

FAO-56 Penman-Monteith reference evapotranspiration (mol H2O m^-2 s^-1), the
numerator of the aridity index `AI = PET_annual/precip_annual` behind `f0`. The
`f0(AI)` relation of Zhou et al. (2025) was fitted with this definition of PET:

    λE = [Δ Rn + ρ_a c_p D / r_a] / [Δ + γ (1 + r_s/r_a)],
    Rn = (1 - α_ref) SW_d + ϵ_sfc (LW_d - σ T^4),

with `Δ` the slope of the saturation vapour pressure curve, `γ` the psychrometric
constant, `D` the vapour pressure deficit, and `r_a = 208/u_2`, `r_s = 70 s m^-1` the
resistances of the 0.12 m reference crop; `u_2` is the wind speed adjusted to 2 m.
`α_ref = 0.23` and `r_s` define that reference surface, so they are not taken from the
simulated canopy.

The ground heat flux is zero (the FAO-56 daily convention; this feeds a yearly total),
and `λE` rather than `Rn` is clipped at zero, so a negative night-time radiative term
offsets the aerodynamic term instead of being dropped.
"""
function potential_evaporation(
    SW_d::FT,
    LW_d::FT,
    T_air::FT,
    P_air::FT,
    q_air::FT,
    u_air::FT,
    h_atmos::FT,
    ϵ_sfc::FT,
    σ::FT,
    M_w::FT,
    thermo_params,
) where {FT}
    α_ref = FT(FAO56_ALBEDO)
    r_s = FT(FAO56_SURFACE_RESISTANCE)
    Rn = (1 - α_ref) * SW_d + ϵ_sfc * (LW_d - σ * T_air^4)

    λv = TP.LH_v0(thermo_params)
    R_v = TP.R_v(thermo_params)
    q = max(q_air, zero(FT))
    c_p = TP.cp_d(thermo_params) * (1 - q) + TP.cp_v(thermo_params) * q

    # Clausius-Clapeyron slope de_sat/dT, and the psychrometric constant with the
    # dry-to-vapour gas constant ratio standing in for the molar mass ratio.
    e_sat = Thermodynamics.saturation_vapor_pressure(
        thermo_params,
        T_air,
        Thermodynamics.Liquid(),
    )
    Δ = e_sat * λv / (R_v * T_air^2)
    γ = c_p * P_air * R_v / (TP.R_d(thermo_params) * λv)

    D = Thermodynamics.vapor_pressure_deficit(
        thermo_params,
        T_air,
        P_air,
        q_air,
    )
    ρ_a = Thermodynamics.air_density(thermo_params, T_air, P_air, q_air)

    # Wind at the reference 2 m (FAO-56 Eq. 47); the relation is anchored on the
    # reference crop, so heights below it are held at 2 m rather than extrapolated.
    u_2 =
        u_air * FT(FAO56_WIND_A) /
        log(FT(FAO56_WIND_B) * max(h_atmos, FT(2)) - FT(FAO56_WIND_C))
    r_a = FT(FAO56_RA_WIND) / max(u_2, sqrt(eps(FT)))

    λE = (Δ * Rn + ρ_a * c_p * D / r_a) / (Δ + γ * (1 + r_s / r_a))
    return max(λE, zero(FT)) / (λv * M_w)
end

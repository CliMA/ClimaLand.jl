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
    """Fraction of annual precipitation available for transpiration (dimensionless, 0-1).
    Following Zhou et al. (2025), f0 = 0.65 at the energy-water limitation transition.
    In arid regions, f0 can be lower: f0 = 0.65 * exp(-0.604 * ln^2(AI/1.9)) where AI is aridity index.
    Default value 0.65 assumes optimal water use efficiency."""
    f0::FT
    """Long-term memory timescale (s) of the A0 and precipitation running-sum totals
    that set LAI_max and the steady-state LAI. Default 1 year; a longer value filters
    the seasonal cycle more strongly."""
    tau_long_term::FT
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
        f0 = FT(toml_dict["optimal_lai_f0"]),
        tau_long_term = FT(toml_dict["optimal_lai_tau_long_term"]),
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


"""
    compute_A0_inst!(dst, p, canopy)

Write the instantaneous potential GPP `A0` (mol CO2 m^-2 s^-1) into `dst`, computed
with fAPAR = 1 from the current drivers and soil-moisture stress. Called once per
cache update; both the daily and annual A0 running sums (`A0_daily`, `A0_annual`)
read the result, so `A0` is evaluated once per step rather than once per variable.

Uses air temperature (not canopy temperature), since `A0` is a potential GPP under
reference conditions, independent of canopy energy-balance feedbacks.
"""
function compute_A0_inst!(dst, p, canopy)
    FT = eltype(canopy.biomass.parameters)
    pmodel_parameters = canopy.photosynthesis.parameters
    pmodel_constants = canopy.photosynthesis.constants
    fractional_c3 = canopy.photosynthesis.fractional_c3
    earth_param_set = canopy.earth_param_set
    # VPD clipped away from zero (the P-model divides by sqrt(VPD)).
    VPD = @. lazy(
        max(
            Thermodynamics.vapor_pressure_deficit(
                LP.thermodynamic_parameters(earth_param_set),
                p.drivers.T,
                p.drivers.P,
                p.drivers.q,
            ),
            sqrt(eps(FT)),
        ),
    )
    PPFD = @. lazy(
        compute_PPFD(
            p.canopy.radiative_transfer.par_d,
            canopy.radiative_transfer.parameters.λ_γ_PAR,
            pmodel_constants.lightspeed,
            pmodel_constants.planck_h,
            pmodel_constants.N_a,
        ),
    )
    @. dst =
        compute_A0_daily(
            fractional_c3,
            pmodel_parameters,
            pmodel_constants,
            p.drivers.T,
            p.drivers.P,
            VPD,
            p.drivers.c_co2,
            PPFD,
            p.canopy.soil_moisture_stress.βm,
        ) / pmodel_constants.Mc
    return nothing
end

"""
    compute_LAI_target!(dst, Y, p, canopy)

Write the instantaneous steady-state LAI target `L_opt` (Zhou et al. 2025 Eqs.
11-15) into `dst`, from the prognostic trailing daily and annual potential-GPP
totals (`A0_daily`, `A0_annual`) and annual precipitation (`precip_annual`) in `Y`,
plus the growing-season water-limitation inputs held in the cache. The prognostic
`LAI` (a `RunningMean`) relaxes toward this target in the biomass
`compute_exp_tendency!`, applying the Eq. 16 acclimation lag continuously through
the time-stepper.
"""
function compute_LAI_target!(dst, Y, p, canopy)
    parameters = canopy.biomass.parameters
    pmodel_parameters = canopy.photosynthesis.parameters
    pmodel_constants = canopy.photosynthesis.constants
    fractional_c3 = canopy.photosynthesis.fractional_c3
    ca = p.drivers.c_co2
    P_air = p.drivers.P
    # chi uses the growing-season mean VPD, not the instantaneous VPD.
    chi = @. lazy(
        compute_chi(
            pmodel_parameters,
            pmodel_constants,
            p.drivers.T,
            P_air,
            p.canopy.biomass.vpd_gs,
            ca,
            fractional_c3,
        ),
    )
    @. dst = compute_L_steady_target(
        Y.canopy.biomass.A0_daily,
        parameters.k,
        Y.canopy.biomass.A0_annual,
        parameters.z,
        p.canopy.biomass.GSL,
        parameters.sigma,
        Y.canopy.biomass.precip_annual,
        p.canopy.biomass.f0,
        ca * P_air,  # ca_pa: CO2 partial pressure (Pa)
        chi,
        p.canopy.biomass.vpd_gs,
    )
    return nothing
end

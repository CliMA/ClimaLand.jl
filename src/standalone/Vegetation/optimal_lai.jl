export AbstractLAIModel,
    OptimalLAIModel,
    OptimalLAIParameters,
    update_LAI!,
    initialize_LAI!,
    set_historical_cache!,
    make_OptimalLAI_callback,
    call_update_optimal_LAI

"""
    AbstractLAIModel{FT <: AbstractFloat}

An abstract type for LAI (Leaf Area Index) model parameterizations.
"""
abstract type AbstractLAIModel{FT <: AbstractFloat} <:
              AbstractCanopyComponent{FT} end

"""
    OptimalLAIParameters{FT<:AbstractFloat}

The required parameters for the optimal LAI model based on Zhou et al. (2025).

Water limitation is handled through the f₀×P/A₀ term following Zhou et al. (2025) Equation 11,
where P is annual precipitation and A₀ is annual potential GPP (computed with β=1).

# References
Zhou et al. (2025) "A General Model for the Seasonal to Decadal Dynamics of Leaf Area"
Global Change Biology. https://onlinelibrary.wiley.com/doi/pdf/10.1111/gcb.70125

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct OptimalLAIParameters{FT <: AbstractFloat}
    """Light extinction coefficient (dimensionless), typically 0.5"""
    k::FT
    """Unit cost of constructing and maintaining leaves (mol m⁻² yr⁻¹), globally fitted as 12.227 mol m⁻² yr⁻¹"""
    z::FT
    """Dimensionless parameter representing departure from square-wave LAI dynamics, globally fitted as 0.771"""
    sigma::FT
    """Smoothing factor for exponential moving average (dimensionless, 0-1). Set to 0.067 for ~15 days of memory"""
    alpha::FT
    """Fraction of annual precipitation available for transpiration (dimensionless, 0-1).
    Following Zhou et al. (2025), f₀ = 0.65 at the energy-water limitation transition.
    In arid regions, f₀ can be lower: f₀ = 0.65 × exp(-0.604 × ln²(AI/1.9)) where AI is aridity index.
    Default value 0.65 assumes optimal water use efficiency."""
    f0::FT
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
    )
end

"""
    OptimalLAIModel{FT, OLPT <: OptimalLAIParameters{FT}, GD} <: AbstractLAIModel{FT}

An implementation of the optimal LAI model from Zhou et al. (2025).

This model predicts seasonal to decadal dynamics of leaf area index based on
optimality principles, balancing energy and water constraints.

# Fields
- `parameters`: Required parameters for the optimal LAI model
- `gsl_a0_data`: NamedTuple with spatially varying GSL (growing season length in days)
  and A0_annual (annual potential GPP in mol CO₂ m⁻² yr⁻¹) fields, typically created
  using `optimal_lai_initial_conditions`.

# References
Zhou et al. (2025) "A General Model for the Seasonal to Decadal Dynamics of Leaf Area"
Global Change Biology. https://onlinelibrary.wiley.com/doi/pdf/10.1111/gcb.70125
"""
struct OptimalLAIModel{FT, OLPT <: OptimalLAIParameters{FT}, GD} <:
       AbstractLAIModel{FT}
    "Required parameters for the optimal LAI model"
    parameters::OLPT
    "Spatially varying GSL and A0_annual data"
    gsl_a0_data::GD
end

Base.eltype(::OptimalLAIModel{FT, OLPT, GD}) where {FT, OLPT, GD} = FT

"""
    OptimalLAIModel{FT}(parameters::OptimalLAIParameters{FT}, gsl_a0_data)

Outer constructor for the OptimalLAIModel struct.

# Arguments
- `parameters`: OptimalLAIParameters for the model
- `gsl_a0_data`: NamedTuple with spatially varying GSL and A0_annual fields,
  typically created using `optimal_lai_initial_conditions`.
"""
function OptimalLAIModel{FT}(
    parameters::OptimalLAIParameters{FT},
    gsl_a0_data,
) where {FT <: AbstractFloat}
    return OptimalLAIModel{FT, typeof(parameters), typeof(gsl_a0_data)}(
        parameters,
        gsl_a0_data,
    )
end

"""
    ClimaLand.auxiliary_vars(model::OptimalLAIModel)
    ClimaLand.auxiliary_types(model::OptimalLAIModel)
    ClimaLand.auxiliary_domain_names(model::OptimalLAIModel)

Defines the auxiliary variables for the OptimalLAIModel:
- `LAI`: leaf area index (m² m⁻²)
- `A0_daily`: daily potential GPP from previous day (mol CO₂ m⁻² day⁻¹), with actual β
- `A0_annual`: annual potential GPP from previous year (mol CO₂ m⁻² yr⁻¹), with β=1 (no moisture stress)
- `A0_daily_acc`: accumulator for current day's potential GPP with actual β (mol CO₂ m⁻² day⁻¹)
- `A0_annual_acc`: accumulator for current year's potential GPP with β=1 (mol CO₂ m⁻² yr⁻¹)
- `A0_annual_daily_acc`: accumulator for current day's potential GPP with β=1 (mol CO₂ m⁻² day⁻¹)
- `days_since_reset`: counter for days since last A0_annual reset (0-365, resets when reaching 365)
- `GSL`: growing season length (days), spatially varying
- `precip_annual`: mean annual precipitation (mol H₂O m⁻² yr⁻¹), for water limitation in LAI_max
- `vpd_gs`: mean VPD during growing season (Pa), for water limitation WUE factor in LAI_max
"""
ClimaLand.auxiliary_vars(model::OptimalLAIModel) =
    (:LAI, :A0_daily, :A0_annual, :A0_daily_acc, :A0_annual_acc, :A0_annual_daily_acc, :days_since_reset, :GSL, :precip_annual, :vpd_gs)
ClimaLand.auxiliary_types(model::OptimalLAIModel{FT}) where {FT} =
    (FT, FT, FT, FT, FT, FT, FT, FT, FT, FT)
ClimaLand.auxiliary_domain_names(::OptimalLAIModel) =
    (:surface, :surface, :surface, :surface, :surface, :surface, :surface, :surface, :surface, :surface)

ClimaLand.name(::AbstractLAIModel) = :lai_model

"""
    initialize_LAI!(p, model::OptimalLAIModel; initial_LAI, initial_A0_daily, initial_A0_annual)

Initialize the LAI and A0 fields to given initial values.

# Arguments
- `p`: auxiliary state
- `model::OptimalLAIModel`: the optimal LAI model
- `initial_LAI`: initial LAI value (m² m⁻², default 1.0)
- `initial_A0_daily`: initial daily potential GPP (mol CO₂ m⁻² day⁻¹, default 0.5)
- `initial_A0_annual`: initial annual potential GPP (mol CO₂ m⁻² yr⁻¹, default 258.0)
"""
function initialize_LAI!(
    p,
    model::OptimalLAIModel;
    initial_LAI = eltype(model)(1.0),
    initial_A0_daily = eltype(model)(0.5),
    initial_A0_annual = eltype(model)(258.0),
)
    FT = eltype(model)
    p.canopy.lai_model.LAI .= initial_LAI
    p.canopy.lai_model.A0_daily .= initial_A0_daily
    p.canopy.lai_model.A0_annual .= initial_A0_annual
    p.canopy.lai_model.A0_daily_acc .= FT(0)
    p.canopy.lai_model.A0_annual_acc .= FT(0)
    p.canopy.lai_model.A0_annual_daily_acc .= FT(0)
    p.canopy.lai_model.days_since_reset .= FT(0)
end

# Conversion factor from mm yr⁻¹ to mol H₂O m⁻² yr⁻¹
# 1 mm = 1 kg m⁻² = 1000 g m⁻²; Molar mass of water = 18.015 g/mol
# So: 1 mm yr⁻¹ = 1000/18.015 mol m⁻² yr⁻¹ ≈ 55.51 mol m⁻² yr⁻¹
const MM_TO_MOL_H2O = 1000.0 / 18.015

"""
    set_historical_cache!(p, Y0, model::OptimalLAIModel, canopy; A0_daily)

The optimal LAI model requires initialization of LAI and A0 values before the simulation.

GSL, A0_annual, precip_annual, vpd_gs, and lai_init are taken from `model.gsl_a0_data`, which
contains spatially varying fields typically created using `optimal_lai_initial_conditions`.

LAI is initialized from MODIS satellite observations (`lai_init`) rather than computing
equilibrium from model equations. This provides a realistic starting point that reduces
spin-up time and matches observed vegetation patterns.

# Arguments
- `A0_daily`: Initial daily potential GPP with actual β (mol CO₂ m⁻² day⁻¹), default 0.5.

# Notes
- precip_annual is converted from mm yr⁻¹ (input) to mol H₂O m⁻² yr⁻¹ (internal units)
  to match Zhou et al. (2025) formulation for water limitation.
- vpd_gs is the mean VPD during growing season (Pa), used in the WUE factor for water limitation.
- lai_init comes from MODIS first timestep, providing spatially-varying initial conditions.
"""
function set_historical_cache!(
    p,
    Y0,
    model::OptimalLAIModel,
    canopy;
    A0_daily = eltype(model.parameters)(0.5),
)
    parameters = model.parameters
    gsl_a0_data = model.gsl_a0_data
    FT = eltype(parameters)

    # Get spatially varying data from gsl_a0_data
    GSL = gsl_a0_data.GSL
    A0_annual = gsl_a0_data.A0_annual
    precip_annual_mm = gsl_a0_data.precip_annual  # in mm yr⁻¹
    vpd_gs = gsl_a0_data.vpd_gs  # in Pa
    lai_init = gsl_a0_data.lai_init  # MODIS first timestep

    L = p.canopy.lai_model.LAI

    # Initialize LAI from MODIS satellite observations
    # This provides realistic spatially-varying initial conditions that reduce spin-up time
    L .= lai_init

    # Initialize A0 variables (supports both scalar and Field inputs via .=)
    p.canopy.lai_model.A0_daily .= A0_daily
    p.canopy.lai_model.A0_annual .= A0_annual
    p.canopy.lai_model.A0_daily_acc .= FT(0)
    p.canopy.lai_model.A0_annual_acc .= FT(0)
    p.canopy.lai_model.A0_annual_daily_acc .= FT(0)
    p.canopy.lai_model.days_since_reset .= FT(0)

    # Store GSL in the cache (spatially varying field)
    p.canopy.lai_model.GSL .= GSL

    # Store precip_annual converted to mol H₂O m⁻² yr⁻¹
    # (input is in mm yr⁻¹, Zhou et al. uses mol m⁻² yr⁻¹)
    @. p.canopy.lai_model.precip_annual = precip_annual_mm * FT(MM_TO_MOL_H2O)

    # Store vpd_gs (already in Pa)
    p.canopy.lai_model.vpd_gs .= vpd_gs
end

"""
    compute_L_max(Ao_annual, k, z, precip_annual, f0, ca_pa, chi, vpd_gs)

Compute seasonal maximum leaf area index (LAI_max) based on annual potential GPP
and water availability, following Zhou et al. (2025) Equation 11.

LAI_max is determined by the minimum of energy-limited and water-limited fAPAR:
- Energy-limited: fAPAR_energy = 1 - z/(k×A₀)
- Water-limited: fAPAR_water = f₀×P/A₀ × (cₐ(1-χ))/(1.6×D)

# Arguments
- `Ao_annual::FT`: Annual total potential GPP (mol CO₂ m⁻² yr⁻¹), computed with β=1.
- `k::FT`: Light extinction coefficient (dimensionless), typically 0.5
- `z::FT`: Unit cost of constructing and maintaining leaves (mol m⁻² yr⁻¹), 12.227
- `precip_annual::FT`: Mean annual precipitation (mol H₂O m⁻² yr⁻¹)
- `f0::FT`: Fraction of precipitation available for transpiration (dimensionless), 0.65
- `ca_pa::FT`: Ambient CO₂ partial pressure (Pa), typically ~40 Pa at 400 ppm
- `chi::FT`: Optimal ratio of intercellular to ambient CO₂ (dimensionless), typically 0.7-0.8
- `vpd_gs::FT`: Mean vapor pressure deficit during growing season (Pa)

# Returns
- `LAI_max::FT`: Seasonal maximum leaf area index (m² m⁻²)

# Notes
Following Zhou et al. (2025) Equation 11:
```
fAPAR_max = min{1 - z/(k×A₀), f₀×P/A₀ × (cₐ(1-χ))/(1.6×D)}
```
The first term is energy-limited (carbon gain vs leaf cost trade-off).
The second term is water-limited (precipitation constrains transpiration, scaled by
intrinsic water use efficiency iWUE = cₐ(1-χ)/(1.6×D)).

The iWUE factor converts water flux to carbon flux:
- cₐ(1-χ): CO₂ drawdown from ambient to intercellular (Pa)
- 1.6×D: VPD adjusted for CO₂/H₂O diffusivity ratio (Pa)

# References
Zhou et al. (2025) Global Change Biology, Equation 11
"""
function compute_L_max(
    Ao_annual::FT,      # mol CO₂ m⁻² yr⁻¹
    k::FT,              # dimensionless
    z::FT,              # mol m⁻² yr⁻¹
    precip_annual::FT,  # mol H₂O m⁻² yr⁻¹
    f0::FT,             # dimensionless
    ca_pa::FT,          # Pa
    chi::FT,            # dimensionless
    vpd_gs::FT,         # Pa
) where {FT}
    # Handle edge case: very small or zero Ao_annual (e.g., polar regions)
    # When Ao_annual ≈ 0, z / (k * Ao_annual) → ∞, causing numerical issues.
    # Use ifelse for GPU compatibility.
    Ao_annual_safe = ifelse(Ao_annual < eps(FT), eps(FT), Ao_annual)

    # Energy-limited fAPAR (Equation 11, first term)
    # Plants optimize leaf area to maximize carbon gain minus construction cost
    fAPAR_energy = FT(1) - z / (k * Ao_annual_safe)

    # Water-limited fAPAR (Equation 11, second term)
    # fAPAR_water = f₀ × P / A₀ × (cₐ(1-χ)) / (1.6×D)
    # The iWUE factor (cₐ(1-χ))/(1.6×D) converts water flux to carbon flux
    # Guard against zero VPD
    vpd_safe = ifelse(vpd_gs < eps(FT), eps(FT), vpd_gs)
    iWUE_factor = (ca_pa * (FT(1) - chi)) / (FT(1.6) * vpd_safe)
    fAPAR_water = f0 * precip_annual / Ao_annual_safe * iWUE_factor

    # fAPAR_max is the minimum of energy and water constraints (Equation 11)
    fAPAR_max = min(fAPAR_energy, fAPAR_water)

    # Ensure fAPAR is in valid range [0, 1]
    fAPAR_max = max(FT(0), min(FT(1), fAPAR_max))

    # Convert fAPAR to LAI using Beer's law (Equation 12)
    # fAPAR = 1 - exp(-k × LAI)  →  LAI = -(1/k) × ln(1 - fAPAR)
    # Guard against fAPAR_max = 1 which would give -log(0) = Inf
    fAPAR_max_safe = min(fAPAR_max, FT(1) - eps(FT))
    LAI_max = -(FT(1) / k) * log(FT(1) - fAPAR_max_safe)

    return LAI_max
end

fAPAR_max_fun(k::FT, LAI_max::FT) where {FT} = FT(1) - exp(-k * LAI_max)

"""
    compute_chi(is_c3, pmodel_parameters, pmodel_constants, T, P_air, VPD, ca)

Compute the optimal ratio of intercellular to ambient CO₂ (χ = ci/ca) using P-model.

# Arguments
- `is_c3`: Photosynthesis mechanism (1 for C3, 0 for C4)
- `pmodel_parameters`: P-model parameters (including β)
- `pmodel_constants`: P-model constants (including Γstar25, Kc25, Ko25, etc.)
- `T::FT`: Temperature (K)
- `P_air::FT`: Atmospheric pressure (Pa)
- `VPD::FT`: Vapor pressure deficit (Pa)
- `ca::FT`: Ambient CO₂ mixing ratio (mol/mol)

# Returns
- `chi::FT`: Optimal ci/ca ratio (dimensionless), typically 0.7-0.85 for C3, ~0.4 for C4
"""
function compute_chi(
    is_c3::AbstractFloat,
    pmodel_parameters,
    pmodel_constants,
    T::FT,
    P_air::FT,
    VPD::FT,
    ca::FT,
) where {FT}
    return is_c3 > 0.5 ? c3_compute_chi(pmodel_parameters, pmodel_constants, T, P_air, VPD, ca) :
                         c4_compute_chi(pmodel_parameters, pmodel_constants, T, P_air, VPD, ca)
end

"""
    c3_compute_chi(pmodel_parameters, pmodel_constants, T, P_air, VPD, ca)

Compute the optimal ratio of intercellular to ambient CO₂ (χ = ci/ca) for C3 plants using P-model.
"""
function c3_compute_chi(
    pmodel_parameters,
    pmodel_constants,
    T::FT,
    P_air::FT,
    VPD::FT,
    ca::FT,
) where {FT}
    (; β) = pmodel_parameters
    (; R, Kc25, Ko25, To, ΔHkc, ΔHko, Drel, ΔHΓstar, Γstar25, oi, ρ_water) = pmodel_constants

    # Convert ca to partial pressure
    ca_pp = ca * P_air

    # Compute P-model intermediates
    Γstar = co2_compensation_pmodel(T, To, P_air, R, ΔHΓstar, Γstar25)
    ηstar = compute_viscosity_ratio(T, To, ρ_water)
    Kmm = compute_Kmm(T, P_air, Kc25, Ko25, ΔHkc, ΔHko, To, R, oi)

    # Compute ξ (sensitivity to dryness)
    ξ = sqrt(β * (Kmm + Γstar) / (Drel * ηstar))

    # Compute ci and chi
    # Guard against zero/negative VPD
    VPD_safe = max(VPD, eps(FT))
    ci = intercellular_co2_pmodel(ξ, ca_pp, Γstar, VPD_safe)
    chi = ci / ca_pp

    # Clamp chi to valid range [0, 1]
    return clamp(chi, FT(0), FT(1))
end

"""
    c4_compute_chi(pmodel_parameters, pmodel_constants, T, P_air, VPD, ca)

Compute the optimal ratio of intercellular to ambient CO₂ (χ = ci/ca) for C4 plants.

For C4 plants, chi is typically around 0.3-0.4 due to the CO2-concentrating mechanism.
This uses an empirical value based on literature (Ehleringer & Cerling).
"""
function c4_compute_chi(
    pmodel_parameters,
    pmodel_constants,
    T::FT,
    P_air::FT,
    VPD::FT,
    ca::FT,
) where {FT}
    # C4 plants typically have ci/ca around 0.4, relatively constant
    # compared to C3 plants which vary from 0.7-0.85
    return FT(0.4)
end

"""
    compute_m(GSL, LAI_max, Ao_annual, sigma)

Compute the parameter m, which represents the ratio of steady-state LAI to steady-state GPP.

This implements Equation 20 from Zhou et al. (2025). The parameter m quantifies the
relationship between LAI and GPP dynamics, representing the extent to which seasonal LAI
dynamics depart from a "square wave" (where maximum LAI would be maintained throughout
the growing season).

# Arguments
- `GSL::FT`: Growing season length (days). Defined as the length of continuous period
  above 0°C longer than 5 days.
- `LAI_max::FT`: Seasonal maximum leaf area index (m² m⁻², dimensionless)
- `Ao_annual::FT`: Annual total potential GPP (mol m⁻² yr⁻¹). This is the integral of daily
  A₀ over the year.
- `sigma::FT`: Dimensionless parameter representing departure from square-wave LAI dynamics.
  Globally fitted as σ = 0.771
- `k::FT`: Light extinction coefficient (dimensionless)

# Returns
- `m::FT`: Parameter relating steady-state LAI to steady-state GPP (dimensionless, units
  work out as: days × m²m⁻² / (mol m⁻² yr⁻¹ × dimensionless) with implicit conversion)

# Example values
```julia
GSL = 180.0          # days (6-month growing season)
LAI_max = 1.6        # m² m⁻²
Ao_annual = 100.0     # mol m⁻² yr⁻¹
sigma = 0.771        # dimensionless
k = 0.5

m = compute_m(GSL, LAI_max, Ao_annual, sigma, k)
# Returns: ~4 (dimensionless, after unit conversion)
```

# References
Zhou et al. (2025) Global Change Biology, Equation 20
"""
function compute_m(
    GSL::FT,        # days
    LAI_max::FT,    # m² m⁻² (dimensionless)
    Ao_annual::FT,     # mol m⁻² yr⁻¹
    sigma::FT,      # dimensionless
    k::FT,
) where {FT}
    # Equation 20: m = (σ × GSL × LAI_max) / (A₀_sum × fAPAR_max)
    fAPAR_max = fAPAR_max_fun(k, LAI_max)

    # Handle edge case: when LAI_max ≈ 0 (very low productivity regions like poles),
    # fAPAR_max = 1 - exp(-k*0) = 0, causing division by zero → NaN.
    # In this case, return m = 0 since there's no vegetation to support.
    # Use ifelse for GPU compatibility (no branching).
    denominator = Ao_annual * fAPAR_max
    m = ifelse(
        denominator < eps(FT),
        zero(FT),
        (sigma * GSL * LAI_max) / denominator,
    )
    return m
end

const MINARG = -inv(Base.MathConstants.e)

"""
    _lambertw0_initial_guess(x::T) where {T<:AbstractFloat}

Provide a robust initial guess for the Lambert W₀ function for use in iterative solvers.

# Arguments
- `x::T`: Input value, should be ≥ -1/e

# Returns
- Initial guess for W₀(x)

# Algorithm
- For x > 1: uses log(x) - log(log(x)) approximation
- For x < -0.32 (near -1/e): uses series expansion for accurate convergence near branch point
- For -0.32 ≤ x ≤ 1: uses max(x, -0.3) as a simple starting point
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
    lambertw0(x::T; maxiter::Int = 16) where {T<:AbstractFloat}

Compute the principal branch (W₀) of the Lambert W function for x ∈ [-1/e, ∞).

This is a GPU-device-friendly implementation using a fixed number of Halley iterations.
The Lambert W function satisfies W(x)·exp(W(x)) = x.

# Arguments
- `x::T`: Input value, must be ≥ -1/e ≈ -0.36788
- `maxiter::Int`: Maximum number of Halley iterations (default: 16)

# Returns
- `W::T`: Lambert W₀(x), the principal branch value, or NaN for invalid inputs

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
@inline function lambertw0(x::T; maxiter::Int = 16) where {T <: AbstractFloat}
    if !(isfinite(x)) || x < T(MINARG)
        return T(NaN)
    end
    w = _lambertw0_initial_guess(x)
    for i in 1:maxiter
        ew = exp(w)
        f = w * ew - x
        # Halley denominator
        # Special case: when w ≈ -1, both numerator and denominator approach 0
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
- `Ao_daily::FT`: Daily potential GPP (mol m⁻² day⁻¹). This is the GPP that would be
  achieved if fAPAR = 1, calculated from LUE × PPFD.
- `m::FT`: Parameter relating steady-state LAI to steady-state GPP (dimensionless), from
  `compute_m()`
- `k::FT`: Light extinction coefficient (dimensionless), typically 0.5
- `LAI_max::FT`: Seasonal maximum LAI constraint (m² m⁻², dimensionless)

# Returns
- `L_steady::FT`: Steady-state leaf area index (m² m⁻², dimensionless). Always ≥ 0.

# Example values
```julia
Ao_daily = 0.4      # mol m⁻² day⁻¹ (summer day in temperate forest)
m = 7.0           # dimensionless (from compute_m)
k = 0.5             # dimensionless
LAI_max = 3.0       # m² m⁻²

L_steady = compute_steady_state_LAI(Ao_daily, m, k, LAI_max)
# Returns: ~1.4 m² m⁻² (close to LAI_max during peak growing season)
```

# Notes
The solution uses the Lambert W₀ function: L_s = min{μ + (1/k)W₀[-kμ exp(-kμ)], LAI_max}
where μ = m × A₀. The result is constrained to be non-negative and below LAI_max.

# References
Zhou et al. (2025) Global Change Biology, Equations 13-15
"""
function compute_steady_state_LAI(
    Ao_daily::FT,  # mol m⁻² day⁻¹
    m::FT,         # dimensionless
    k::FT,         # dimensionless
    LAI_max::FT,   # m² m⁻² (dimensionless)
) where {FT}
    # μ = m × A₀ (Equation 15)
    mu = m * Ao_daily

    # Compute argument for Lambert W function
    arg = -k * mu * exp(-k * mu)

    # Check if argument is in valid range for W₀ branch: [-1/e, 0]
    # If outside this range, use boundary solution
    if arg < -FT(1) / FT(exp(1))
        # Beyond valid range; use maximum possible LAI
        L_s = LAI_max
    else
        # Equation 15: L_s = μ + (1/k) × W₀[-k μ exp(-k μ)]
        # Using our custom lambertw0 function (W₀ is the principal branch)
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
    update_LAI!(LAI_prev, L_steady, alpha, local_noon_mask)

Update LAI using exponential weighted moving average to represent lag between carbon
allocation and steady-state LAI.

This implements Equation 16 from Zhou et al. (2025). The exponential moving average
represents the time lag (days to months) for photosynthate allocation to leaves and
leaf development.

# Arguments
- `LAI_prev::FT`: LAI from previous time step (m² m⁻², dimensionless)
- `L_steady::FT`: Current steady-state LAI (m² m⁻², dimensionless), from
  `compute_steady_state_LAI()`
- `alpha::FT`: Smoothing factor (dimensionless, 0-1). Set to 0.067 for ~15 days of memory,
  meaning LAI[t] = 0.067 × L_steady[t] + 0.933 × LAI[t-1]
- `local_noon_mask::FT`: A mask (0 or 1) indicating whether the current time is within the local noon window.

# Returns
- `LAI_new::FT`: Updated leaf area index (m² m⁻², dimensionless). Always ≥ 0.

# Example values
```julia
LAI_prev = 3.5      # m² m⁻² (LAI from yesterday)
L_steady = 4.0      # m² m⁻² (steady-state LAI for today)
alpha = 0.067       # dimensionless (15-day memory)
local_noon_mask = 1.0  # update at local noon

LAI_new = update_LAI!(LAI_prev, L_steady, alpha, local_noon_mask)
# Returns: ~3.53 m² m⁻² (slowly approaching steady-state)
```

# Notes
The parameter α controls the response time:
- α = 0.067 ≈ 15 days of memory (paper default)
- α = 0.1 ≈ 10 days of memory (faster response)
- α = 0.033 ≈ 30 days of memory (slower response)

The time scale τ ≈ 1/α days.

# References
Zhou et al. (2025) Global Change Biology, Equation 16
"""
function update_LAI!(
    LAI_prev::FT,   # m² m⁻² (dimensionless)
    L_steady::FT,   # m² m⁻² (dimensionless)
    alpha::FT,      # dimensionless (0-1)
    local_noon_mask::FT,
) where {FT}
    if local_noon_mask == FT(1.0)
        # Equation 16: LAI_sim = α × L_s[t] + (1-α) × LAI_sim[t-1]
        LAI_new = alpha * L_steady + (1 - alpha) * LAI_prev

        # Ensure non-negative (important for numerical stability, especially at initialization)
        LAI_new = max(zero(FT), LAI_new)
        return LAI_new
    else
        return LAI_prev
    end
end

"""
    update_optimal_LAI(local_noon_mask, A0_daily, L, k, A0_annual, z, GSL, sigma, alpha, precip_annual, f0, ca_pa, chi, vpd_gs)

Update LAI using the optimal LAI model with precomputed daily and annual potential GPP.

# Arguments
- `local_noon_mask::FT`: Mask (0 or 1) indicating if it's local noon
- `A0_daily::FT`: Daily potential GPP (mol CO₂ m⁻² day⁻¹), with actual β
- `L::FT`: Current LAI (m² m⁻²)
- `k::FT`: Light extinction coefficient
- `A0_annual::FT`: Annual potential GPP (mol CO₂ m⁻² yr⁻¹), with β=1 (no moisture stress)
- `z::FT`: Unit cost of constructing and maintaining leaves (mol m⁻² yr⁻¹)
- `GSL::FT`: Growing season length (days)
- `sigma::FT`: Dimensionless parameter for LAI dynamics
- `alpha::FT`: Smoothing factor for exponential moving average (~15-day memory)
- `precip_annual::FT`: Mean annual precipitation (mol H₂O m⁻² yr⁻¹)
- `f0::FT`: Fraction of precipitation available for transpiration (dimensionless)
- `ca_pa::FT`: Ambient CO₂ partial pressure (Pa)
- `chi::FT`: Optimal ratio of intercellular to ambient CO₂ (dimensionless)
- `vpd_gs::FT`: Mean vapor pressure deficit during growing season (Pa)

# Returns
Updated LAI value.

# Notes
Following Zhou et al. (2025):
- A0_daily uses actual β (soil moisture stress) to drive daily LAI dynamics
- A0_annual uses β=1 (no moisture stress) for LAI_max computation
- Water limitation enters LAI_max through the f₀×P/A₀ × (cₐ(1-χ))/(1.6×D) term (Equation 11)
"""
function update_optimal_LAI(
    local_noon_mask::FT,
    A0_daily::FT,
    L::FT, # m2 m-2
    k::FT,
    A0_annual::FT, # mol CO2 m-2 y-1, computed with β=1
    z::FT, # mol m-2 yr-1, leaf construction cost
    GSL::FT, # days, growing season length
    sigma::FT, # dimensionless
    alpha::FT, # dimensionless (~15-day memory)
    precip_annual::FT, # mol H2O m-2 yr-1, mean annual precipitation
    f0::FT, # dimensionless, fraction of precip for transpiration
    ca_pa::FT, # Pa, ambient CO2 partial pressure
    chi::FT, # dimensionless, optimal ci/ca ratio
    vpd_gs::FT, # Pa, mean VPD during growing season
) where {FT}
    LAI_max = compute_L_max(A0_annual, k, z, precip_annual, f0, ca_pa, chi, vpd_gs)
    m = compute_m(GSL, LAI_max, A0_annual, sigma, k)
    L_steady = compute_steady_state_LAI(A0_daily, m, k, LAI_max)
    L = update_LAI!(L, L_steady, alpha, local_noon_mask)
    return L
end

"""
    compute_PPFD(par_d::FT, λ_γ_PAR::FT, lightspeed::FT, planck_h::FT, N_a::FT) where {FT}

Compute photosynthetic photon flux density (PPFD) from PAR downwelling flux.
This is the total incoming PAR (not absorbed), in units of mol photons m⁻² s⁻¹.
"""
function compute_PPFD(
    par_d::FT,
    λ_γ_PAR::FT,
    lightspeed::FT,
    planck_h::FT,
    N_a::FT,
) where {FT}
    energy_per_mole_photon_par = planck_h * lightspeed * N_a / λ_γ_PAR
    return par_d / energy_per_mole_photon_par
end

"""
    call_update_optimal_LAI(p, Y, t, current_date; canopy, dt, local_noon)

Updates LAI and accumulates A0 at each timestep. At local noon, finalizes daily A0
and updates LAI. On year change (Jan 1), finalizes annual A0.

The function:
1. Computes instantaneous A0_daily using the P-model with fAPAR=1 and actual β (soil moisture stress)
2. Computes instantaneous A0_annual using the P-model with fAPAR=1 and β=1 (no moisture stress)
3. Accumulates A0_daily and A0_annual to their respective accumulators
4. At local noon: finalizes daily A0, adds A0_annual to annual accumulator, updates LAI
5. On Jan 1: finalizes annual A0

A0_daily uses actual β to drive daily LAI dynamics. A0_annual uses β=1 following Zhou et al.
(2025); water limitation enters LAI_max through the f₀×P/A₀ term using annual precipitation
and VPD (to be implemented).

GSL (Growing Season Length) is read from p.canopy.lai_model.GSL, which supports spatially
varying values initialized via set_historical_cache! or optimal_lai_initial_conditions.
"""
function call_update_optimal_LAI(
    p,
    Y,
    t,
    current_date;
    canopy,
    dt,
    local_noon,
)
    FT = eltype(canopy.lai_model.parameters)

    # Compute local noon mask
    local_noon_mask = @. lazy(get_local_noon_mask(t, dt, local_noon))

    # Get P-model parameters and constants for computing A0
    pmodel_parameters = canopy.photosynthesis.parameters
    pmodel_constants = canopy.photosynthesis.constants
    is_c3 = canopy.photosynthesis.is_c3

    # Get drivers for A0 computation
    T_canopy = canopy_temperature(canopy.energy, canopy, Y, p)
    P_air = p.drivers.P
    ca = p.drivers.c_co2  # mol/mol
    earth_param_set = canopy.earth_param_set

    # Get soil moisture stress factor (β)
    βm = p.canopy.soil_moisture_stress.βm

    # Compute VPD (clipped to avoid numerical issues)
    VPD = @. lazy(
        max(
            ClimaLand.vapor_pressure_deficit(
                p.drivers.T,
                p.drivers.P,
                p.drivers.q,
                LP.thermodynamic_parameters(earth_param_set),
            ),
            sqrt(eps(FT)),
        ),
    )

    # Compute PPFD from PAR downwelling (total incoming, not absorbed)
    par_d = p.canopy.radiative_transfer.par_d
    λ_γ_PAR = canopy.radiative_transfer.parameters.λ_γ_PAR
    PPFD = @. lazy(
        compute_PPFD(
            par_d,
            λ_γ_PAR,
            pmodel_constants.lightspeed,
            pmodel_constants.planck_h,
            pmodel_constants.N_a,
        ),
    )

    # Compute instantaneous A0_daily with actual β (kg C m⁻² s⁻¹)
    # This captures soil moisture stress for daily GPP driving LAI dynamics
    A0_daily_inst = @. lazy(
        compute_A0_daily(
            is_c3,
            pmodel_parameters,
            pmodel_constants,
            T_canopy,
            P_air,
            VPD,
            ca,
            PPFD,
            βm,
        ),
    )

    # Compute instantaneous A0_annual with β=1 (no moisture stress) (kg C m⁻² s⁻¹)
    # Following Zhou et al. (2025), water limitation enters LAI_max through
    # the f₀×P/A₀ term (using annual precipitation and VPD), not through β.
    A0_annual_inst = @. lazy(
        compute_A0_potential(
            is_c3,
            pmodel_parameters,
            pmodel_constants,
            T_canopy,
            P_air,
            VPD,
            ca,
            PPFD,
        ),
    )

    # Convert A0 from kg C m⁻² s⁻¹ to mol CO₂ m⁻² per timestep
    # Mc = 0.0120107 kg/mol C
    dt_seconds = FT(float(dt))
    Mc = pmodel_constants.Mc
    A0_daily_per_timestep = @. lazy(A0_daily_inst * dt_seconds / Mc)
    A0_annual_per_timestep = @. lazy(A0_annual_inst * dt_seconds / Mc)

    # Accumulate to daily A0 (with actual β) - drives L_steady
    @. p.canopy.lai_model.A0_daily_acc += A0_daily_per_timestep

    # Accumulate to daily A0_annual (with β=1) - will be added to annual accumulator at noon
    @. p.canopy.lai_model.A0_annual_daily_acc += A0_annual_per_timestep

    # Get parameters from the LAI model
    parameters = canopy.lai_model.parameters

    # At local noon: finalize daily A0, update annual accumulator, update LAI
    # Inline the update logic to avoid tuple broadcasting issues.
    # Compute finalized A0 values before any resets.
    # Year reset happens when days_since_reset >= 365 (ensures full year of accumulation)
    A0_daily_final = @. p.canopy.lai_model.A0_daily_acc
    A0_annual_daily_final = @. p.canopy.lai_model.A0_annual_daily_acc
    A0_annual_final = @. ifelse(
        local_noon_mask == FT(1) &&
        p.canopy.lai_model.days_since_reset >= FT(365),
        p.canopy.lai_model.A0_annual_acc,
        p.canopy.lai_model.A0_annual,
    )

    # Compute ca_pa (CO₂ partial pressure in Pa) and chi for water limitation term
    # chi is computed using P-model formula with growing season mean VPD (vpd_gs)
    ca_pa = @. lazy(ca * P_air)  # Convert mol/mol to Pa
    chi = @. lazy(
        compute_chi(
            is_c3,
            pmodel_parameters,
            pmodel_constants,
            T_canopy,
            P_air,
            p.canopy.lai_model.vpd_gs,  # Use growing season mean VPD
            ca,
        ),
    )

    # Update LAI at noon using the finalized values
    # GSL, precip_annual, and vpd_gs are read from the cache (supports spatially varying values)
    @. p.canopy.lai_model.LAI = ifelse(
        local_noon_mask == FT(1),
        update_optimal_LAI(
            FT(1),  # At noon, we always update
            A0_daily_final,
            p.canopy.lai_model.LAI,
            parameters.k,
            A0_annual_final, # already accounts for year change
            parameters.z,
            p.canopy.lai_model.GSL,
            parameters.sigma,
            parameters.alpha,
            p.canopy.lai_model.precip_annual,
            parameters.f0,
            ca_pa,
            chi,
            p.canopy.lai_model.vpd_gs,
        ),
        p.canopy.lai_model.LAI,
    )

    # Update A0_annual if year reset (365 days accumulated)
    @. p.canopy.lai_model.A0_annual = ifelse(
        local_noon_mask == FT(1) &&
        p.canopy.lai_model.days_since_reset >= FT(365),
        A0_annual_final,
        p.canopy.lai_model.A0_annual,
    )

    # Update A0_daily at noon (this stores the daily value with actual β)
    @. p.canopy.lai_model.A0_daily = ifelse(
        local_noon_mask == FT(1),
        A0_daily_final,
        p.canopy.lai_model.A0_daily,
    )

    # Reset A0_annual_acc if year reset, otherwise add new daily A0_annual (β=1) value at noon
    # On reset: start fresh with today's value (A0_annual_daily_final)
    @. p.canopy.lai_model.A0_annual_acc = ifelse(
        local_noon_mask == FT(1) &&
        p.canopy.lai_model.days_since_reset >= FT(365),
        A0_annual_daily_final,  # Reset and add today's value
        ifelse(
            local_noon_mask == FT(1),
            p.canopy.lai_model.A0_annual_acc + A0_annual_daily_final,
            p.canopy.lai_model.A0_annual_acc,
        ),
    )

    # Reset A0_daily_acc at noon
    @. p.canopy.lai_model.A0_daily_acc = ifelse(
        local_noon_mask == FT(1),
        FT(0),
        p.canopy.lai_model.A0_daily_acc,
    )

    # Reset A0_annual_daily_acc at noon
    @. p.canopy.lai_model.A0_annual_daily_acc = ifelse(
        local_noon_mask == FT(1),
        FT(0),
        p.canopy.lai_model.A0_annual_daily_acc,
    )

    # Update days_since_reset at noon
    # On year reset: start counting from 1 (today is day 1 of new accumulation period)
    # Otherwise: increment by 1
    @. p.canopy.lai_model.days_since_reset = ifelse(
        local_noon_mask == FT(1) &&
        p.canopy.lai_model.days_since_reset >= FT(365),
        FT(1),  # Reset counter, today is day 1
        ifelse(
            local_noon_mask == FT(1),
            p.canopy.lai_model.days_since_reset + FT(1),
            p.canopy.lai_model.days_since_reset,
        ),
    )
end

"""
    update_A0_and_LAI_at_noon(local_noon_mask, A0_daily, A0_annual, A0_daily_acc, A0_annual_acc, A0_annual_daily_acc, days_since_reset, L, k, z, GSL, sigma, alpha, precip_annual, f0, ca_pa, chi, vpd_gs)

At local noon: finalize daily A0, add to annual, reset daily accumulator, and update LAI.
On year reset (days_since_reset >= 365): finalize annual A0.

Note: A0_daily_acc contains daily GPP with actual β (soil moisture stress).
      A0_annual_daily_acc contains daily GPP with β=1 (for LAI_max computation).

Returns tuple: (A0_daily, A0_annual, A0_daily_acc, A0_annual_acc, A0_annual_daily_acc, days_since_reset, L)
"""
function update_A0_and_LAI_at_noon(
    local_noon_mask::FT,
    A0_daily::FT,
    A0_annual::FT,
    A0_daily_acc::FT,
    A0_annual_acc::FT,
    A0_annual_daily_acc::FT,
    days_since_reset::FT,
    L::FT,
    k::FT,
    z::FT,
    GSL::FT,
    sigma::FT,
    alpha::FT,
    precip_annual::FT,
    f0::FT,
    ca_pa::FT,
    chi::FT,
    vpd_gs::FT,
) where {FT}
    if local_noon_mask == FT(1)
        # Check for year reset (365 days accumulated = full year)
        year_reset = days_since_reset >= FT(365)

        if year_reset
            # Finalize annual A0 before starting new accumulation period
            A0_annual = A0_annual_acc
        end

        # Finalize daily A0 (with actual β for daily LAI dynamics)
        A0_daily = A0_daily_acc
        A0_daily_acc = FT(0)

        # Add today's A0_annual (with β=1, for LAI_max) to annual accumulator
        # On year reset, start fresh with today's value
        if year_reset
            A0_annual_acc = A0_annual_daily_acc
            days_since_reset = FT(1)  # Today is day 1 of new period
        else
            A0_annual_acc += A0_annual_daily_acc
            days_since_reset += FT(1)
        end
        A0_annual_daily_acc = FT(0)

        # Update LAI using the optimal LAI model
        L = update_optimal_LAI(
            local_noon_mask,
            A0_daily,
            L,
            k,
            A0_annual,
            z,
            GSL,
            sigma,
            alpha,
            precip_annual,
            f0,
            ca_pa,
            chi,
            vpd_gs,
        )

        return A0_daily, A0_annual, A0_daily_acc, A0_annual_acc, A0_annual_daily_acc, days_since_reset, L
    else
        return A0_daily, A0_annual, A0_daily_acc, A0_annual_acc, A0_annual_daily_acc, days_since_reset, L
    end
end

"""
    make_OptimalLAI_callback(::Type{FT}, t0::ITime, dt, canopy; longitude) where {FT <: AbstractFloat}

This constructs an IntervalBasedCallback for the optimal LAI model that:
1. Computes and accumulates potential GPP (A₀) at each timestep
2. Updates LAI using an exponential moving average at local noon
3. Tracks daily and annual A₀ sums

We check for local noon using the provided `longitude` every dt.
The time of local noon is expressed in seconds UTC and neglects the effects of obliquity and eccentricity, so
it is constant throughout the year.

# Arguments
- `FT`: The floating-point type used in the model (e.g., `Float32`, `Float64`).
- `t0`: ITime, with epoch in UTC.
- `dt`: timestep
- `canopy`: the canopy object containing the optimal LAI model parameters.
- `longitude`: optional longitude in degrees for local noon calculation (default is `nothing`, which means
    that it will be inferred from the canopy domain).

# Notes
- Daily A₀ is computed with fAPAR=1 and actual β (soil moisture stress) - drives L_steady
- Annual A₀ is computed with fAPAR=1 and actual β (soil moisture stress) - used for LAI_max
- Water limitation is captured through β in both daily and annual A₀, allowing vegetation
  structure to adapt on annual timescales while responding to climate change and flushing events
- Daily A₀ is accumulated over each day and finalized at local noon
- Annual A₀ is accumulated and reset on January 1
- GSL (Growing Season Length) is read from p.canopy.lai_model.GSL, which supports spatially
  varying values initialized via set_historical_cache! or optimal_lai_initial_conditions.
"""
function make_OptimalLAI_callback(
    ::Type{FT},
    t0,
    dt,
    canopy;
    longitude = nothing,
) where {FT <: AbstractFloat}
    function seconds_after_midnight(d)
        return FT(
            Hour(d).value * 3600 +
            Minute(d).value * 60 +
            Second(d).value,
        )
    end

    if isnothing(longitude)
        try
            longitude = get_long(canopy.domain.space.surface)
        catch e
            error(
                "Longitude must be provided explicitly if the domain you are working on does not \
                  have axes that specify longitude $e",
            )
        end
    end

    # this computes the time of local noon in seconds UTC without considering the
    # effects of obliquity and orbital eccentricity, so it is constant throughout the year
    # the max error is on the order of 20 minutes
    seconds_in_a_day = IP.day(IP.InsolationParameters(FT))
    start_t = seconds_after_midnight(date(t0))
    start_date = date(t0)
    local_noon = @. seconds_in_a_day * (FT(1 / 2) - longitude / 360) # allocates, but only on init

    affect! =
        (integrator) -> begin
            # Compute current date from t0 and elapsed time
            elapsed_days = floor(Int, float(integrator.t) / seconds_in_a_day)
            current_date = start_date + Dates.Day(elapsed_days)

            call_update_optimal_LAI(
                integrator.p,
                integrator.u,
                (float(integrator.t) + start_t) % (seconds_in_a_day), # current time in seconds UTC
                current_date;
                canopy = canopy,
                dt = dt,
                local_noon = local_noon,
            )
        end

    return IntervalBasedCallback(
        dt,         # period of this callback
        t0,         # simulation start
        dt,         # integration timestep
        affect!;
    )
end

"""
    get_model_callbacks(component::OptimalLAIModel, canopy; t0, Δt)

Creates the optimal LAI callback and returns it as a single element tuple of model callbacks.

# Notes
- Daily A₀ is computed with actual β (soil moisture stress) from the soil moisture stress model
- Annual A₀ is computed with actual β (soil moisture stress) for LAI_max computation
- Water limitation enters through β in both daily and annual A₀, allowing vegetation structure
  to adapt on annual timescales while responding to climate change and flushing events
- GSL (Growing Season Length) is read from p.canopy.lai_model.GSL, which supports spatially
  varying values initialized via set_historical_cache! or optimal_lai_initial_conditions.
"""
function get_model_callbacks(
    component::OptimalLAIModel{FT},
    canopy;
    t0,
    Δt,
) where {FT}
    lai_cb = make_OptimalLAI_callback(
        FT,
        t0,
        Δt,
        canopy,
    )
    return (lai_cb,)
end

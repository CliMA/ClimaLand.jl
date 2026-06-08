export update_albedo!,
    CLMTwoBandSoilAlbedo, ConstantTwoBandSoilAlbedo, CompositionBasedSoilAlbedo

abstract type AbstractSoilAlbedoParameterization end

"""
    ConstantTwoBandSoilAlbedo{
        SF <: Union{AbstractFloat, ClimaCore.Fields.Field},
    } <: AbstractSoilAlbedoParameterization

A parameterization for soil albedo: the soil albedo is
defined in two bands (PAR and NIR), can spatially vary or
be set to scalar, but is assumed not to vary in time (or
with water content).
"""
struct ConstantTwoBandSoilAlbedo{
    SF <: Union{AbstractFloat, ClimaCore.Fields.Field},
} <: AbstractSoilAlbedoParameterization
    "Soil PAR Albedo"
    PAR_albedo::SF
    "Soil NIR Albedo"
    NIR_albedo::SF
end

function ConstantTwoBandSoilAlbedo{FT}(;
    PAR_albedo::FT = FT(0.2),
    NIR_albedo::FT = FT(0.4),
) where {FT}
    ConstantTwoBandSoilAlbedo{FT}(PAR_albedo, NIR_albedo)
end


"""
    update_albedo!(bc::AtmosDrivenFluxBC, albedo::ConstantTwoBandSoilAlbedo, p, soil_domain, model_parameters)

Updates PAR and NIR albedo using the temporally-constant parameters
provided in `albedo`; these values may be spatially varying.

For the temporally-constant albedo model, there is no need to update
the values each step, or have an allocated spot in the cache for them. This
can be optimized in the future.
"""
function update_albedo!(
    bc::AtmosDrivenFluxBC,
    albedo::ConstantTwoBandSoilAlbedo,
    p,
    soil_domain,
    model_parameters,
)
    @. p.soil.PAR_albedo = albedo.PAR_albedo
    @. p.soil.NIR_albedo = albedo.NIR_albedo
end

"""
     CLMTwoBandSoilAlbedo{
        FT <: AbstractFloat,
        SF <: Union{FT, ClimaCore.Fields.Field},
    } <: AbstractSoilAlbedoParameterization

A parameterization for soil albedo: the soil albedo is
defined in two bands (PAR and NIR), and can spatially vary or
be set to scalar. However, it varies temporally due to a
dependence on soil water content at the surface,
via the effective saturation S(θ\\_sfc):
α = α\\_wet*S + α\\_dry*(1-S)

We use a value for θ\\_sfc averaged over the depth `albedo_calc_top_thickness`.
If the model resolution is such that the first layer is thicker than this depth,
the value from the first layer is used.

CLM reference: Lawrence, P.J., and Chase, T.N. 2007. Representing a MODIS consistent land surface in the Community Land Model
(CLM 3.0). J. Geophys. Res. 112:G01023. DOI:10.1029/2006JG000168.
"""
struct CLMTwoBandSoilAlbedo{
    FT <: AbstractFloat,
    SF <: Union{FT, ClimaCore.Fields.Field},
} <: AbstractSoilAlbedoParameterization
    "Soil PAR Albedo dry"
    PAR_albedo_dry::SF
    "Soil NIR Albedo dry"
    NIR_albedo_dry::SF
    "Soil PAR Albedo wet"
    PAR_albedo_wet::SF
    "Soil NIR Albedo wet"
    NIR_albedo_wet::SF
    "Thickness of top of soil used in albedo calculations (m)"
    albedo_calc_top_thickness::FT
end

function CLMTwoBandSoilAlbedo{FT}(;
    PAR_albedo_dry,
    NIR_albedo_dry,
    PAR_albedo_wet,
    NIR_albedo_wet,
    albedo_calc_top_thickness = FT(0.02),
) where {FT}
    return CLMTwoBandSoilAlbedo{FT, typeof(PAR_albedo_dry)}(
        PAR_albedo_dry,
        NIR_albedo_dry,
        PAR_albedo_wet,
        NIR_albedo_wet,
        albedo_calc_top_thickness,
    )
end


"""
    albedo_from_moisture(S_sfc::FT, θ_int::FT, albedo_wet::FT)

Calculates pointwise albedo for any band as a function of soil surface moisture given
the dry and wet albedo values for that band using the CLM parameterization.

CLM reference: Lawrence, P.J., and Chase, T.N. 2007. Representing a MODIS consistent land surface in the Community Land Model
(CLM 3.0). J. Geophys. Res. 112:G01023. DOI:10.1029/2006JG000168.
"""
function albedo_from_moisture(
    S_sfc::FT,
    albedo_dry::FT,
    albedo_wet::FT,
) where {FT}
    return (1 - S_sfc) * albedo_dry + S_sfc * albedo_wet
end


"""
    update_albedo!(bc::AtmosDrivenFluxBC, albedo::CLMTwoBandSoilAlbedo, p, soil_domain, model_parameters)

Calculates and updates PAR and NIR albedo as a function of volumetric soil water content at
the top of the soil. If the soil layers are larger than the specified `albedo_calc_top_thickness`,
the water content of the top layer is used in the calclulation. For the PAR and NIR bands,

α\\_band = α\\_{band,dry} * (1 - S\\_e) +  α\\_{band,wet} * (S\\_e)

where S\\_e is the relative soil wetness above some depth, `albedo_calc_top_thickness`. This
is a modified version of Equation (1) of:

Braghiere, R. K., Wang, Y., Gagné-Landmann, A., Brodrick, P. G., Bloom, A. A., Norton,
A. J., et al. (2023). The importance of hyperspectral soil albedo information for improving
Earth system model projections. AGU Advances, 4, e2023AV000910. https://doi.org/10.1029/2023AV000910

where effective saturation is used in place of volumetric soil water content. The dry and wet
albedo values come from a global soil color map and soil color to albedo map from CLM.

CLM reference: Lawrence, P.J., and Chase, T.N. 2007. Representing a MODIS consistent land surface in the Community Land Model
(CLM 3.0). J. Geophys. Res. 112:G01023. DOI:10.1029/2006JG000168.
"""
function update_albedo!(
    bc::AtmosDrivenFluxBC,
    albedo::CLMTwoBandSoilAlbedo,
    p,
    soil_domain,
    model_parameters,
)
    (;
        PAR_albedo_dry,
        NIR_albedo_dry,
        PAR_albedo_wet,
        NIR_albedo_wet,
        albedo_calc_top_thickness,
    ) = albedo
    FT = eltype(soil_domain.fields.Δz_top)
    # checks if there is at least 1 layer centered within the top soil depth
    if soil_domain.fields.Δz_min < albedo_calc_top_thickness
        # We compute ∫H_θ_dz / N, where N =∫H_dz
        # is a normalization,
        # and then compute θ_sfc = ∫[(H θ)/N]dz.
        # N is integral of 1 from (surface-albedo_calc_top_thickness) to surface
        N = p.soil.sfc_scratch
        # zero all centers lower than boundary, set everything above to one
        @. p.soil.sub_sfc_scratch = ClimaLand.heaviside(
            albedo_calc_top_thickness + sqrt(eps(FT)),
            soil_domain.fields.z_sfc - soil_domain.fields.z,
        )
        ClimaCore.Operators.column_integral_definite!(N, p.soil.sub_sfc_scratch)
        # zeros all effective saturation at levels centered lower than boundary
        @. p.soil.sub_sfc_scratch =
            ClimaLand.heaviside(
                albedo_calc_top_thickness + sqrt(eps(FT)),
                soil_domain.fields.z_sfc - soil_domain.fields.z,
            ) * p.soil.θ_l / N
        # Now use scratch space to compute normed integral
        θ_sfc = p.soil.sfc_scratch
        ClimaCore.Operators.column_integral_definite!(
            θ_sfc,
            p.soil.sub_sfc_scratch,
        )
    else
        # in the case where no layer is centered above boundary, use the values of the top layer
        θ_sfc = ClimaLand.Domains.top_center_to_surface(p.soil.θ_l)
    end
    ν_sfc = ClimaLand.Domains.top_center_to_surface(model_parameters.ν)
    θ_r_sfc = ClimaLand.Domains.top_center_to_surface(model_parameters.θ_r)
    S_sfc = @. lazy(effective_saturation(ν_sfc, θ_sfc, θ_r_sfc))
    @. p.soil.PAR_albedo =
        albedo_from_moisture(S_sfc, PAR_albedo_dry, PAR_albedo_wet)
    @. p.soil.NIR_albedo =
        albedo_from_moisture(S_sfc, NIR_albedo_dry, NIR_albedo_wet)
end


"""
    CompositionBasedSoilAlbedo{FT <: AbstractFloat} <: AbstractSoilAlbedoParameterization

A data-driven parameterization for soil albedo that computes dry albedo from
soil composition and texture using logistic regression, then applies nonlinear
moisture darkening.

This model uses:
- Organic matter (ν_ss_om): darkening effect
- Van Genuchten n parameter: pore size distribution, related to texture
- Coarse fragments (ν_ss_gravel): gravel/rock content

(TODO: possibly add more variables)

The dry albedo uses logistic regression:
    η = η₀ + c_om * ν_ss_om + c_vgn * vg_n + c_cf * ν_ss_gravel
    α_dry = α_min + (α_max - α_min) * σ(η)

where σ(η) = 1/(1 + exp(-η)) is the logistic function that maps (-∞, ∞) → (0, 1).
This ensures α_dry ∈ (α_min, α_max) with smooth, continuous gradients everywhere.

The moisture dependence uses a nonlinear (power-law) relationship:
    α = α_dry * (1 - f_wet * S_e^β)

where S_e is effective saturation, β controls the nonlinearity (β < 1 gives
rapid initial darkening as observed), and f_wet is the maximum moisture
darkening factor.

This model achieves **RMSE = 0.074** and **R² = 0.51** on desert regions,
a 30% RMSE reduction compared to the spatial-mean albedo (TODO: Update numbers). 
The improvement comes from:
- vg_n correlates strongly with albedo (r = 0.71) as it captures texture
- Coarse fragments (r = 0.35) indicate rocky/gravelly surfaces

# Calibration
Default coefficients fitted on bareground albedo in geographic desert regions 
(Sahara, Arabian, Australian, Gobi, Kalahari) using SoilGrids composition and van 
Genuchten soil parameters.

# Physical basis
- **Organic matter**: Strongly absorbs visible light (dark color), reduces albedo
- **Van Genuchten n**: Higher n indicates sandier soils (brighter); lower n
  indicates clay-like soils (darker)
- **Coarse fragments**: Rocky/gravelly surfaces tend to be brighter
- **Moisture**: Water in pore spaces increases absorption; effect is strongest
  at low moisture levels (Lobell & Asner, 2002)

# References
- Lobell, D. B., & Asner, G. P. (2002). Moisture effects on soil reflectance.
  Soil Science Society of America Journal, 66(3), 722-727.

  # TODO: Make all of these coefficients calibratable ClimaParams
"""
struct CompositionBasedSoilAlbedo{FT <: AbstractFloat} <:
       AbstractSoilAlbedoParameterization
    "Base log-odds for PAR"
    η₀_PAR::FT
    "Base log-odds for NIR"
    η₀_NIR::FT
    "Organic matter coefficient for PAR (negative, darkening)"
    c_om_PAR::FT
    "Organic matter coefficient for NIR"
    c_om_NIR::FT
    "Van Genuchten n coefficient for PAR (positive, coarser = brighter)"
    c_vgn_PAR::FT
    "Van Genuchten n coefficient for NIR"
    c_vgn_NIR::FT
    "Coarse fragments coefficient for PAR (positive, rocky = brighter)"
    c_cf_PAR::FT
    "Coarse fragments coefficient for NIR"
    c_cf_NIR::FT
    "Maximum moisture darkening factor"
    f_wet::FT
    "Moisture sensitivity exponent"
    β::FT
    "Minimum allowed albedo"
    α_min::FT
    "Maximum allowed albedo"
    α_max::FT
    "Thickness of top of soil used in albedo calculations (m)"
    albedo_calc_top_thickness::FT
end

"""
    CompositionBasedSoilAlbedo{FT}(;
        η₀_PAR = FT(-2.25),
        η₀_NIR = FT(-2.3),
        c_om_PAR = FT(-0.13),
        c_om_NIR = FT(-0.14),
        c_vgn_PAR = FT(1.24),
        c_vgn_NIR = FT(1.28),
        c_cf_PAR = FT(0.15),
        c_cf_NIR = FT(0.16),
        f_wet = FT(0.50),
        β = FT(0.5),
        α_min = FT(0.04),
        α_max = FT(0.60),
        albedo_calc_top_thickness = FT(0.02),
    ) where {FT}

Construct a `CompositionBasedSoilAlbedo` with calibrated coefficients.

# Calibration details
Coefficients fitted via logistic regression on:
- Target: Bareground shortwave albedo (TODO: this is based on CERES; update with MODIS)
- Predictors: ν_ss_om (SoilGrids), vg_n (van Genuchten), ν_ss_gravel (SoilGrids)
- Training: Geographic desert regions (Sahara, Arabian, Australian, Gobi, Kalahari)
- Validation RMSE: 0.074 (30% reduction vs climatology)
- Validation R²: 0.51

# Expected albedo values

**High vg_n desert (sandy)**: vg_n ≈ 2.5, cf ≈ 0.3, om ≈ 0
- η ≈ -2.25 + 1.24*2.5 + 0.15*0.3 = 0.90 → α ≈ 0.44

**Rocky desert**: vg_n ≈ 1.5, cf ≈ 0.6, om ≈ 0
- η ≈ -2.25 + 1.24*1.5 + 0.15*0.6 = -0.29 → α ≈ 0.28

**Organic soil**: vg_n ≈ 1.5, cf ≈ 0.1, om ≈ 0.1
- η ≈ -2.25 + 1.24*1.5 - 0.13*0.1 + 0.15*0.1 = -0.37 → α ≈ 0.27

# Arguments
- `η₀_PAR`, `η₀_NIR`: Base log-odds (intercept in logistic model)
- `c_om_PAR`, `c_om_NIR`: Organic matter coefficients (negative = darkening)
- `c_vgn_PAR`, `c_vgn_NIR`: Van Genuchten n coefficients (positive = sandier = brighter)
- `c_cf_PAR`, `c_cf_NIR`: Coarse fragments coefficients (positive = brighter)
- `f_wet`: Maximum fractional reduction in albedo when fully saturated
- `β`: Moisture sensitivity exponent (β < 1 for nonlinear rapid initial darkening)
- `α_min`, `α_max`: Bounds of sigmoid output (physical albedo limits)
- `albedo_calc_top_thickness`: Depth for surface moisture averaging (m)
"""
function CompositionBasedSoilAlbedo{FT}(;
    # TODO: All of these coefficients should become calibratable ClimaParams
    η₀_PAR = FT(-2.25),
    η₀_NIR = FT(-2.3),
    c_om_PAR = FT(-0.13),
    c_om_NIR = FT(-0.14),
    c_vgn_PAR = FT(1.24),
    c_vgn_NIR = FT(1.28),
    c_cf_PAR = FT(0.15),
    c_cf_NIR = FT(0.16),
    f_wet = FT(0.50),
    β = FT(0.5),
    α_min = FT(0.04),
    α_max = FT(0.60),
    albedo_calc_top_thickness = FT(0.02),
) where {FT}
    return CompositionBasedSoilAlbedo{FT}(
        η₀_PAR,
        η₀_NIR,
        c_om_PAR,
        c_om_NIR,
        c_vgn_PAR,
        c_vgn_NIR,
        c_cf_PAR,
        c_cf_NIR,
        f_wet,
        β,
        α_min,
        α_max,
        albedo_calc_top_thickness,
    )
end


"""
    dry_albedo_from_composition(
        η₀::FT, c_om::FT, c_vgn::FT, c_cf::FT,
        ν_ss_om::FT, vg_n::FT, ν_ss_gravel::FT,
        α_min::FT, α_max::FT
    ) where {FT}

Compute dry soil albedo from soil composition using logistic regression.

η = η₀ + c_om * ν_ss_om + c_vgn * vg_n + c_cf * ν_ss_gravel
α_dry = α_min + (α_max - α_min) / (1 + exp(-η))

This ensures α_dry ∈ (α_min, α_max) with smooth, continuous derivatives.
"""
function dry_albedo_from_composition(
    η₀::FT,
    c_om::FT,
    c_vgn::FT,
    c_cf::FT,
    ν_ss_om::FT,
    vg_n::FT,
    ν_ss_gravel::FT,
    α_min::FT,
    α_max::FT,
) where {FT}
    η = η₀ + c_om * ν_ss_om + c_vgn * vg_n + c_cf * ν_ss_gravel
    # Numerically stable sigmoid
    σ = if η >= FT(0)
        FT(1) / (FT(1) + exp(-η))
    else
        exp_η = exp(η)
        exp_η / (FT(1) + exp_η)
    end
    return α_min + (α_max - α_min) * σ
end


"""
    albedo_with_nonlinear_moisture(α_dry::FT, S_e::FT, f_wet::FT, β::FT) where {FT}

Apply nonlinear moisture darkening to dry albedo.

α = α_dry * (1 - f_wet * S_e^β)

where β < 1 gives rapid initial darkening (most moisture effect at low saturation),
matching observations that soil darkens quickly with initial wetting.
"""
function albedo_with_nonlinear_moisture(
    α_dry::FT,
    S_e::FT,
    f_wet::FT,
    β::FT,
) where {FT}
    # Clamp S_e to [0, 1] for safety
    S_e_clamped = clamp(S_e, FT(0), FT(1))
    moisture_factor = S_e_clamped^β
    return α_dry * (FT(1) - f_wet * moisture_factor)
end


"""
    update_albedo!(bc::AtmosDrivenFluxBC, albedo::CompositionBasedSoilAlbedo, p, soil_domain, model_parameters)

Computes and updates PAR and NIR albedo based on soil composition and moisture.

The dry albedo is computed via logistic regression on soil composition:
    η = η₀ + c_om * ν_ss_om + c_vgn * vg_n + c_cf * ν_ss_gravel
    α_dry = α_min + (α_max - α_min) * σ(η)

where σ(η) = 1/(1 + exp(-η)) is the logistic function.

Then moisture darkening is applied nonlinearly:
    α = α_dry * (1 - f_wet * S_e^β)

This approach links albedo directly to the soil properties available in the model.
"""
function update_albedo!(
    bc::AtmosDrivenFluxBC,
    albedo::CompositionBasedSoilAlbedo,
    p,
    soil_domain,
    model_parameters,
)
    (;
        η₀_PAR,
        η₀_NIR,
        c_om_PAR,
        c_om_NIR,
        c_vgn_PAR,
        c_vgn_NIR,
        c_cf_PAR,
        c_cf_NIR,
        f_wet,
        β,
        α_min,
        α_max,
        albedo_calc_top_thickness,
    ) = albedo
    FT = eltype(soil_domain.fields.Δz_top)

    # Get surface soil properties
    ν_ss_om_sfc =
        ClimaLand.Domains.top_center_to_surface(model_parameters.ν_ss_om)
    ν_ss_gravel_sfc =
        ClimaLand.Domains.top_center_to_surface(model_parameters.ν_ss_gravel)
    hydrology_cm_sfc =
        ClimaLand.Domains.top_center_to_surface(model_parameters.hydrology_cm)

    # Extract vg_n from the van Genuchten closure
    vg_n_sfc = @. lazy(hydrology_cm_sfc.n)

    # Compute dry albedo from composition
    α_dry_PAR = @. lazy(
        dry_albedo_from_composition(
            η₀_PAR,
            c_om_PAR,
            c_vgn_PAR,
            c_cf_PAR,
            ν_ss_om_sfc,
            vg_n_sfc,
            ν_ss_gravel_sfc,
            α_min,
            α_max,
        ),
    )
    α_dry_NIR = @. lazy(
        dry_albedo_from_composition(
            η₀_NIR,
            c_om_NIR,
            c_vgn_NIR,
            c_cf_NIR,
            ν_ss_om_sfc,
            vg_n_sfc,
            ν_ss_gravel_sfc,
            α_min,
            α_max,
        ),
    )

    # Compute surface soil moisture
    if soil_domain.fields.Δz_min < albedo_calc_top_thickness
        N = p.soil.sfc_scratch
        @. p.soil.sub_sfc_scratch = ClimaLand.heaviside(
            albedo_calc_top_thickness + sqrt(eps(FT)),
            soil_domain.fields.z_sfc - soil_domain.fields.z,
        )
        ClimaCore.Operators.column_integral_definite!(N, p.soil.sub_sfc_scratch)
        @. p.soil.sub_sfc_scratch =
            ClimaLand.heaviside(
                albedo_calc_top_thickness + sqrt(eps(FT)),
                soil_domain.fields.z_sfc - soil_domain.fields.z,
            ) * p.soil.θ_l / N
        θ_sfc = p.soil.sfc_scratch
        ClimaCore.Operators.column_integral_definite!(
            θ_sfc,
            p.soil.sub_sfc_scratch,
        )
    else
        θ_sfc = ClimaLand.Domains.top_center_to_surface(p.soil.θ_l)
    end

    ν_sfc = ClimaLand.Domains.top_center_to_surface(model_parameters.ν)
    θ_r_sfc = ClimaLand.Domains.top_center_to_surface(model_parameters.θ_r)
    S_sfc = @. lazy(effective_saturation(ν_sfc, θ_sfc, θ_r_sfc))

    # Apply nonlinear moisture darkening
    @. p.soil.PAR_albedo =
        albedo_with_nonlinear_moisture(α_dry_PAR, S_sfc, f_wet, β)
    @. p.soil.NIR_albedo =
        albedo_with_nonlinear_moisture(α_dry_NIR, S_sfc, f_wet, β)
end


"""
    update_albedo!(bc::AbstractEnergyHydrologyBC, _...)

Does nothing for boundary conditions where albedo is not used.
"""
function update_albedo!(bc::AbstractEnergyHydrologyBC, _...) end


"""
    ClimaLand.surface_albedo(
        model::EnergyHydrology{FT},
        Y,
        p,
    ) where {FT}

Returns the surface albedo field of the
`EnergyHydrology` soil model.
"""
function ClimaLand.surface_albedo(model::EnergyHydrology{FT}, Y, p) where {FT}
    return @. lazy((p.soil.PAR_albedo + p.soil.NIR_albedo) / 2)
end

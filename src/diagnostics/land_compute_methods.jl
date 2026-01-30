"""
    @diagnostic_compute name model compute

Macro generating a function to compute a land diagnostic,
needed in the ClimaDiagnostics framework.

See the ClimaDiagnostics docmentation for more information.

For example,
@diagnostic_compute "soil_net_radiation" SoilCanopyModel p.soil.R_n

generates the function

function compute_soil_net_radiation!(out, Y, p, t, land_model::SoilCanopyModel)
    if isnothing(out)
        out = zeros(axes(p.soil.R_n)) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        out .= p.soil.R_n # set the land values only since this type of broadcasting respects the mask
        return out
    else
        out .= p.soil.R_n
    end
end

Please note that if a land/sea mask is employed, the values
over the ocean are set to NaN.
"""
macro diagnostic_compute(name, model, compute)
    function_name = Symbol("compute_", name, "!")
    return esc(
        quote
            @with_error function $function_name(
                out,
                Y,
                p,
                t,
                land_model::$model,
            )
                if isnothing(out)
                    out = zeros(axes($compute)) # Allocates
                    fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
                    out .= $compute # set the land values only since this type of broadcasting respects the mask
                    return out
                else
                    out .= $compute
                end
            end
        end,
    )
end

## Helper functions so that we can use the same methods for integrated
## and standalone models

get_canopy(m::Union{SoilCanopyModel, LandModel}) = m.canopy
get_canopy(m::CanopyModel) = m

get_soil(m::Union{SoilCanopyModel, LandModel, SoilSnowModel}) = m.soil
get_soil(m::EnergyHydrology) = m

### Conservation ##
@diagnostic_compute "water_volume_per_area" EnergyHydrology p.soil.total_water
@diagnostic_compute "energy_per_area" EnergyHydrology p.soil.total_energy
@diagnostic_compute "water_volume_per_area_change" EnergyHydrology Y.soil.∫F_vol_liq_water_dt
@diagnostic_compute "energy_per_area_change" EnergyHydrology Y.soil.∫F_e_dt


### BucketModel ###

# variables stored in p (diagnostics variables stored in the cache)
@diagnostic_compute "sw_albedo" BucketModel p.bucket.α_sfc
@diagnostic_compute "latent_heat_flux" BucketModel p.bucket.turbulent_fluxes.lhf
@diagnostic_compute "net_radiation" BucketModel p.bucket.R_n
@diagnostic_compute "sensible_heat_flux" BucketModel p.bucket.turbulent_fluxes.shf
@diagnostic_compute "surface_air_density" BucketModel p.bucket.ρ_sfc
@diagnostic_compute "specific_humidity" BucketModel p.bucket.q_sfc
@diagnostic_compute "surface_temperature" BucketModel p.bucket.T_sfc
@diagnostic_compute "vapor_flux" BucketModel p.bucket.turbulent_fluxes.vapor_flux

# variables stored in Y (prognostic or state variables)
@diagnostic_compute "snow_water_equivalent" BucketModel Y.bucket.σS
@diagnostic_compute "soil_temperature" BucketModel Y.bucket.T
@diagnostic_compute "subsurface_water_storage" BucketModel Y.bucket.W
@diagnostic_compute "surface_water_content" BucketModel Y.bucket.Ws

### Union{SoilCanopyModel, LandModel} ###

# variables stored in p (diagnostics variables stored in the cache)

## Canopy Module ##

# Canopy - Solar Induced Fluorescence
@diagnostic_compute "solar_induced_fluorescence" Union{
    SoilCanopyModel,
    LandModel,
    CanopyModel,
} p.canopy.sif.SIF

# Canopy - Autotrophic respiration
@diagnostic_compute "autotrophic_respiration" Union{
    SoilCanopyModel,
    LandModel,
    CanopyModel,
} p.canopy.autotrophic_respiration.Ra

# Net Ecosystem Exchange (NEE = ER - GPP)
function compute_net_ecosystem_exchange!(
    out,
    Y,
    p,
    t,
    land_model::Union{SoilCanopyModel{FT}, LandModel{FT}},
) where {FT}
    # Compute ER first
    if isnothing(out)
        out = zeros(land_model.soil.domain.space.surface)
        fill!(field_values(out), NaN)
    end
    compute_total_respiration!(out, Y, p, t, land_model)
    @. out -= p.canopy.photosynthesis.GPP
    return out
end

# Canopy - Conductance
function compute_stomatal_conductance!(
    out,
    Y,
    p,
    t,
    land_model::Union{SoilCanopyModel, LandModel, CanopyModel},
)
    canopy = get_canopy(land_model)
    conductance_model = canopy.conductance
    compute_stomatal_conductance!(out, Y, p, t, canopy, conductance_model)
end


function compute_stomatal_conductance!(
    out,
    Y,
    p,
    t,
    canopy,
    conductance_model::MedlynConductanceModel,
)
    (; g1, g0, Drel) = conductance_model.parameters
    earth_param_set = canopy.earth_param_set
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    An_leaf = get_An_leaf(p, canopy.photosynthesis)
    if isnothing(out)
        out = zeros(canopy.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        @. out = medlyn_conductance(
            g0,
            Drel,
            medlyn_term.(
                g1,
                p.drivers.T,
                p.drivers.P,
                p.drivers.q,
                thermo_params,
            ),
            An_leaf,
            p.drivers.c_co2,
        )
        return out
    else
        @. out = medlyn_conductance(
            g0,
            Drel,
            medlyn_term(
                g1,
                p.drivers.T,
                p.drivers.P,
                p.drivers.q,
                thermo_params,
            ),
            An_leaf,
            p.drivers.c_co2,
        )
    end
end

function compute_stomatal_conductance!(
    out,
    Y,
    p,
    t,
    canopy,
    conductance_model::PModelConductance,
)
    (; Drel) = conductance_model.parameters
    c_co2_air = p.drivers.c_co2
    P_air = p.drivers.P
    ci = p.canopy.photosynthesis.ci             # internal CO2 partial pressure, Pa
    An_leaf = get_An_leaf(p, canopy.photosynthesis)
    if isnothing(out)
        out = zeros(canopy.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        @. out =
            gs_h2o_pmodel(ci / (c_co2_air * P_air), c_co2_air, An_leaf, Drel)
        return out
    else
        @. out =
            gs_h2o_pmodel(ci / (c_co2_air * P_air), c_co2_air, An_leaf, Drel)
    end
end

function compute_canopy_transpiration!(
    out,
    Y,
    p,
    t,
    land_model::Union{CanopyModel, SoilCanopyModel, LandModel},
)
    canopy = get_canopy(land_model)
    # Convert to a mass flux by multiplying by the density of liquid
    # water
    if isnothing(out)
        out = zeros(canopy.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        @. out = p.canopy.turbulent_fluxes.vapor_flux * 1000
        return out
    else
        @. out = p.canopy.turbulent_fluxes.vapor_flux * 1000
    end
end

# Canopy - Energy
@diagnostic_compute "canopy_latent_heat_flux" Union{
    SoilCanopyModel,
    LandModel,
    CanopyModel,
} p.canopy.turbulent_fluxes.lhf
@diagnostic_compute "canopy_sensible_heat_flux" Union{
    SoilCanopyModel,
    LandModel,
    CanopyModel,
} p.canopy.turbulent_fluxes.shf
@diagnostic_compute "latent_heat_flux" CanopyModel p.canopy.turbulent_fluxes.lhf
@diagnostic_compute "sensible_heat_flux" CanopyModel p.canopy.turbulent_fluxes.shf

# Canopy - Hydraulics
function compute_leaf_water_potential!(
    out,
    Y,
    p,
    t,
    land_model::Union{CanopyModel, SoilCanopyModel, LandModel},
)
    canopy = get_canopy(land_model)
    hydraulics = canopy.hydraulics
    n_stem = hydraulics.n_stem
    n_leaf = hydraulics.n_leaf
    n = n_stem + n_leaf
    if isnothing(out)
        out = zeros(canopy.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        out .= p.canopy.hydraulics.ψ.:($n)
        return out
    else
        out .= p.canopy.hydraulics.ψ.:($n)
    end
end

# @diagnostic_compute "flux_per_ground_area" Union{SoilCanopyModel, LandModel} p.canopy.hydraulics.fa # return a Tuple
@diagnostic_compute "root_flux_per_ground_area" Union{
    SoilCanopyModel,
    LandModel,
    CanopyModel,
} p.canopy.hydraulics.fa_roots
@diagnostic_compute "leaf_area_index" Union{
    SoilCanopyModel,
    LandModel,
    CanopyModel,
} p.canopy.biomass.area_index.leaf

# Canopy - Soil moisture stress
@diagnostic_compute "moisture_stress_factor" Union{
    SoilCanopyModel,
    LandModel,
    CanopyModel,
} p.canopy.soil_moisture_stress.βm

# Canopy - Hydraulics
@diagnostic_compute "root_area_index" Union{
    SoilCanopyModel,
    LandModel,
    CanopyModel,
} p.canopy.biomass.area_index.root
@diagnostic_compute "stem_area_index" Union{
    SoilCanopyModel,
    LandModel,
    CanopyModel,
} p.canopy.biomass.area_index.stem

# Canopy - Photosynthesis
@diagnostic_compute "photosynthesis_gross_canopy" Union{
    SoilCanopyModel,
    LandModel,
    CanopyModel,
} p.canopy.photosynthesis.GPP

@diagnostic_compute "photosynthesis_net_leaf" Union{SoilCanopyModel, LandModel} get_An_leaf(
    p,
    land_model.canopy.photosynthesis,
)
@diagnostic_compute "photosynthesis_net_leaf" CanopyModel get_An_leaf(
    p,
    land_model.photosynthesis,
)
@diagnostic_compute "respiration_leaf" Union{SoilCanopyModel, LandModel} get_Rd_leaf(
    p,
    land_model.canopy.photosynthesis,
)
@diagnostic_compute "respiration_leaf" CanopyModel get_Rd_leaf(
    p,
    land_model.photosynthesis,
)
@diagnostic_compute "vcmax25" Union{SoilCanopyModel, LandModel} get_Vcmax25_leaf(
    p,
    land_model.canopy.photosynthesis,
)
@diagnostic_compute "vcmax25" CanopyModel get_Vcmax25_leaf(
    p,
    land_model.photosynthesis,
)

# Canopy - Radiative Transfer
@diagnostic_compute "near_infrared_radiation_down" Union{
    SoilCanopyModel,
    LandModel,
    CanopyModel,
} p.canopy.radiative_transfer.nir_d
@diagnostic_compute "near_infrared_radiation_absorbed" Union{
    SoilCanopyModel,
    LandModel,
    CanopyModel,
} p.canopy.radiative_transfer.nir.abs
@diagnostic_compute "near_infrared_radiation_reflected" Union{
    SoilCanopyModel,
    LandModel,
    CanopyModel,
} p.canopy.radiative_transfer.nir.refl
@diagnostic_compute "near_infrared_radiation_transmitted" Union{
    SoilCanopyModel,
    LandModel,
    CanopyModel,
} p.canopy.radiative_transfer.nir.trans
@diagnostic_compute "photosynthetically_active_radiation_down" Union{
    SoilCanopyModel,
    LandModel,
    CanopyModel,
} p.canopy.radiative_transfer.par_d
@diagnostic_compute "photosynthetically_active_radiation_absorbed" Union{
    SoilCanopyModel,
    LandModel,
    CanopyModel,
} p.canopy.radiative_transfer.par.abs
@diagnostic_compute "photosynthetically_active_radiation_reflected" Union{
    SoilCanopyModel,
    LandModel,
    CanopyModel,
} p.canopy.radiative_transfer.par.refl
@diagnostic_compute "photosynthetically_active_radiation_transmitted" Union{
    SoilCanopyModel,
    LandModel,
    CanopyModel,
} p.canopy.radiative_transfer.par.trans
@diagnostic_compute "radiation_longwave_net" Union{
    SoilCanopyModel,
    LandModel,
    CanopyModel,
} p.canopy.radiative_transfer.LW_n
@diagnostic_compute "radiation_shortwave_net" Union{
    SoilCanopyModel,
    LandModel,
    CanopyModel,
} p.canopy.radiative_transfer.SW_n

## Drivers Module ##

@diagnostic_compute "soil_organic_carbon" Union{
    SoilCanopyModel,
    LandModel,
    SoilCO2Model,
} Y.soilco2.SOC # SOC is now prognostic
@diagnostic_compute "pressure" Union{SoilCanopyModel, LandModel, CanopyModel} p.drivers.P
@diagnostic_compute "rainfall" Union{SoilCanopyModel, LandModel, CanopyModel} p.drivers.P_liq
@diagnostic_compute "radiation_longwave_down" Union{
    SoilCanopyModel,
    LandModel,
    CanopyModel,
} p.drivers.LW_d
@diagnostic_compute "radiation_shortwave_down" Union{
    SoilCanopyModel,
    LandModel,
    CanopyModel,
} p.drivers.SW_d
@diagnostic_compute "snowfall" Union{SoilCanopyModel, LandModel, CanopyModel} p.drivers.P_snow
@diagnostic_compute "tair" Union{SoilCanopyModel, LandModel, CanopyModel} p.drivers.T
@diagnostic_compute "specific_humidity" Union{
    SoilCanopyModel,
    LandModel,
    CanopyModel,
} p.drivers.q
@diagnostic_compute "wind_speed" Union{SoilCanopyModel, LandModel, CanopyModel} p.drivers.u

function compute_precip!(
    out,
    Y,
    p,
    t,
    land_model::Union{EnergyHydrology, CanopyModel, SoilCanopyModel, LandModel},
)
    soil = get_soil(land_model)
    if isnothing(out)
        out = zeros(soil.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        @. out = (p.drivers.P_liq + p.drivers.P_snow) * 1000 # density of liquid water (1000kg/m^3)
        return out
    else
        @. out = (p.drivers.P_liq + p.drivers.P_snow) * 1000# density of liquid water (1000kg/m^3)
    end
end

## Soil Module ##

@diagnostic_compute "infiltration" Union{
    EnergyHydrology,
    SoilCanopyModel,
    LandModel,
} p.soil.infiltration
@diagnostic_compute "soil_hydraulic_conductivity" Union{
    SoilCanopyModel,
    LandModel,
    EnergyHydrology,
} p.soil.K
@diagnostic_compute "soil_thermal_conductivity" Union{
    SoilCanopyModel,
    LandModel,
    EnergyHydrology,
} p.soil.κ
@diagnostic_compute "soil_water_potential" Union{
    SoilCanopyModel,
    LandModel,
    EnergyHydrology,
} p.soil.ψ
@diagnostic_compute "soil_temperature" Union{
    SoilCanopyModel,
    LandModel,
    EnergyHydrology,
} p.soil.T
@diagnostic_compute "soil_net_radiation" Union{
    SoilCanopyModel,
    LandModel,
    EnergyHydrology,
} p.soil.R_n

function compute_10cm_water_mass!(
    out,
    Y,
    p,
    t,
    land_model::Union{EnergyHydrology{FT}, SoilCanopyModel{FT}, LandModel{FT}},
) where {FT}
    soil = get_soil(land_model)
    ∫Hθdz = p.soil.sfc_scratch
    Hθ = p.soil.sub_sfc_scratch
    z = land_model.soil.domain.fields.z
    depth = FT(-0.1)
    earth_param_set = soil.parameters.earth_param_set
    _ρ_liq = LP.ρ_cloud_liq(earth_param_set)
    _ρ_ice = LP.ρ_cloud_ice(earth_param_set)
    # Convert from volumetric water content to water mass per unit volume using density
    @. Hθ = (p.soil.θ_l * _ρ_liq + Y.soil.θ_i * _ρ_ice) * heaviside(z, depth)
    column_integral_definite!(∫Hθdz, Hθ)

    # The layering of the soil model may not coincide with 10 cm exactly, and this could lead
    # to the integral above not exactly representing 10cm.
    # To adjust, divide by the ∫heaviside(z, depth) dz, and then multiply by 10cm
    H = p.subsfc_scratch
    @. H = heaviside(z, depth)
    ∫Hdz = p.sfc_scratch
    column_integral_definite!(∫Hdz, H)

    if isnothing(out)
        out = zeros(soil.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        @. out = ∫Hθdz / ∫Hdz * FT(0.1)
        return out
    else
        @. out = ∫Hθdz / ∫Hdz * FT(0.1)
    end
end
function compute_soil_albedo!(
    out,
    Y,
    p,
    t,
    land_model::Union{SoilCanopyModel{FT}, LandModel{FT}, EnergyHydrology{FT}},
) where {FT}
    soil = get_soil(land_model)
    if isnothing(out)
        out = zeros(soil.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        @. out = (p.soil.PAR_albedo + p.soil.NIR_albedo) / 2
        return out
    else
        @. out = (p.soil.PAR_albedo + p.soil.NIR_albedo) / 2
    end
end

# Soil - Turbulent Fluxes
@diagnostic_compute "soil_latent_heat_flux" Union{
    SoilCanopyModel,
    LandModel,
    EnergyHydrology,
} p.soil.turbulent_fluxes.lhf
@diagnostic_compute "soil_sensible_heat_flux" Union{
    SoilCanopyModel,
    LandModel,
    EnergyHydrology,
} p.soil.turbulent_fluxes.shf
@diagnostic_compute "vapor_flux" Union{
    SoilCanopyModel,
    LandModel,
    EnergyHydrology,
} p.soil.turbulent_fluxes.vapor_flux_liq # should add ice here

# Soil - SoilCO2
function compute_heterotrophic_respiration!(
    out,
    Y,
    p,
    t,
    land_model::SoilCO2Model{FT},
) where {FT}
    if isnothing(out)
        out = zeros(land_model.soil.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        @. out = p.top_bc * FT(83.26)
        return out
    else
        out .= p.top_bc .* FT(83.26)
    end
end # Convert from kg C to mol CO2.
function compute_heterotrophic_respiration!(
    out,
    Y,
    p,
    t,
    land_model::Union{SoilCanopyModel{FT}, LandModel{FT}},
) where {FT}
    if isnothing(out)
        out = zeros(land_model.soil.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        @. out = p.soilco2.top_bc * FT(83.26)
        return out
    else
        out .= p.soilco2.top_bc .* FT(83.26)
    end
end # Convert from kg C to mol CO2.
# To convert from kg C to mol CO2, we need to multiply by:
# [3.664 kg CO2/ kg C] x [10^3 g CO2/ kg CO2] x [1 mol CO2/44.009 g CO2] = 83.26 mol CO2/kg C

@diagnostic_compute "soilco2_diffusivity" SoilCO2Model p.soilco2.D
@diagnostic_compute "soilco2_source_microbe" SoilCO2Model p.soilco2.Sm
@diagnostic_compute "soilo2_diffusivity" SoilCO2Model p.soilco2.D_o2
@diagnostic_compute "soilco2_diffusivity" Union{SoilCanopyModel, LandModel} p.soilco2.D
@diagnostic_compute "soilco2_source_microbe" Union{SoilCanopyModel, LandModel} p.soilco2.Sm
@diagnostic_compute "soilo2_diffusivity" Union{SoilCanopyModel, LandModel} p.soilco2.D_o2

## Other ##
@diagnostic_compute "sw_albedo" Union{SoilCanopyModel, LandModel} p.α_sfc
@diagnostic_compute "lw_up" Union{SoilCanopyModel, LandModel} p.LW_u
@diagnostic_compute "sw_up" Union{SoilCanopyModel, LandModel} p.SW_u
function compute_sw_up!(out, Y, p, t, land_model::BucketModel)
    α_sfc = ClimaLand.surface_albedo(land_model, Y, p)
    SW_d = p.drivers.SW_d

    if isnothing(out)
        out = zeros(land_model.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        @. out = α_sfc * SW_d
        return out
    else
        @. out = α_sfc * SW_d
    end
end

function compute_lw_up!(out, Y, p, t, land_model::BucketModel)
    LW_d = p.drivers.LW_d
    earth_param_set = land_model.parameters.earth_param_set
    _σ = LP.Stefan(earth_param_set)
    T_sfc = ClimaLand.component_temperature(land_model, Y, p)
    ϵ_sfc = ClimaLand.surface_emissivity(land_model, Y, p)
    if isnothing(out)
        out = zeros(land_model.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        @. out = (1 - ϵ_sfc) * LW_d + ϵ_sfc * _σ * T_sfc^4
        return out
    else
        @. out = (1 - ϵ_sfc) * LW_d + ϵ_sfc * _σ * T_sfc^4
    end
end

function compute_soil_fsat!(
    out,
    Y,
    p,
    t,
    land_model::Union{EnergyHydrology, SoilCanopyModel, LandModel},
)
    soil = get_soil(land_model)
    if isnothing(out)
        out = zeros(soil.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
    end
    out .= Runoff.get_soil_fsat(
        soil.boundary_conditions.top.runoff,
        Y,
        p,
        soil.domain.fields.depth,
    )
    return out
end

function compute_surface_runoff!(
    out,
    Y,
    p,
    t,
    land_model::Union{EnergyHydrology, SoilCanopyModel, LandModel},
)
    soil = get_soil(land_model)
    if isnothing(out)
        out = zeros(soil.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
    end
    out .= Runoff.get_surface_runoff(soil.boundary_conditions.top.runoff, Y, p)
    return out
end

function compute_subsurface_runoff!(
    out,
    Y,
    p,
    t,
    land_model::Union{EnergyHydrology, SoilCanopyModel, LandModel},
)
    soil = get_soil(land_model)
    if isnothing(out)
        out = zeros(soil.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
    end
    out .=
        Runoff.get_subsurface_runoff(soil.boundary_conditions.top.runoff, Y, p)
    return out
end

function compute_saturated_height!(
    out,
    Y,
    p,
    t,
    land_model::Union{EnergyHydrology, SoilCanopyModel, LandModel},
)
    soil = get_soil(land_model)
    if isnothing(out)
        out = zeros(soil.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
    end

    out .=
        Runoff.get_saturated_height(soil.boundary_conditions.top.runoff, Y, p)
    return out
end

function compute_infiltration_capacity!(
    out,
    Y,
    p,
    t,
    land_model::Union{EnergyHydrology, SoilCanopyModel, LandModel},
)
    soil = get_soil(land_model)
    if isnothing(out)
        out = zeros(soil.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
    end
    out .= Runoff.soil_infiltration_capacity(soil, Y, p)
    return out
end

@diagnostic_compute "bottom_water_flux" Union{
    EnergyHydrology,
    SoilCanopyModel,
    LandModel,
} p.soil.bottom_bc.water

function compute_evapotranspiration!(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel{FT},
) where {FT}
    if isnothing(out)
        out = zeros(land_model.soil.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        @. out =
            (
                p.soil.turbulent_fluxes.vapor_flux_liq +
                p.soil.turbulent_fluxes.vapor_flux_ice +
                p.canopy.turbulent_fluxes.vapor_flux
            ) * 1000 # density of liquid water (1000kg/m^3)
        return out
    else
        out .=
            (
                p.soil.turbulent_fluxes.vapor_flux_liq .+
                p.soil.turbulent_fluxes.vapor_flux_ice .+
                p.canopy.turbulent_fluxes.vapor_flux
            ) .* 1000 # density of liquid water (1000kg/m^3)
    end
end

function compute_evapotranspiration!(
    out,
    Y,
    p,
    t,
    land_model::LandModel{FT},
) where {FT}
    if isnothing(out)
        out = zeros(land_model.canopy.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        @. out =
            (
                (1 - p.snow.snow_cover_fraction) *
                p.soil.turbulent_fluxes.vapor_flux_liq +
                (1 - p.snow.snow_cover_fraction) *
                p.soil.turbulent_fluxes.vapor_flux_ice +
                p.canopy.turbulent_fluxes.vapor_flux +
                p.snow.snow_cover_fraction * p.snow.turbulent_fluxes.vapor_flux
            ) * 1000 # density of liquid water (1000kg/m^3)
        return out
    else
        @. out =
            (
                (1 - p.snow.snow_cover_fraction) *
                p.soil.turbulent_fluxes.vapor_flux_liq +
                (1 - p.snow.snow_cover_fraction) *
                p.soil.turbulent_fluxes.vapor_flux_ice +
                p.canopy.turbulent_fluxes.vapor_flux +
                p.snow.snow_cover_fraction * p.snow.turbulent_fluxes.vapor_flux
            ) * 1000 # density of liquid water (1000kg/m^3)
    end
end
@diagnostic_compute "evapotranspiration" CanopyModel p.canopy.turbulent_fluxes.vapor_flux

function compute_total_respiration!(
    out,
    Y,
    p,
    t,
    land_model::Union{SoilCanopyModel{FT}, LandModel{FT}},
) where {FT}
    if isnothing(out)
        out = zeros(land_model.canopy.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        if isnothing(land_model.soilco2)
            @. out = p.canopy.autotrophic_respiration.Ra
        else
            @. out =
                p.soilco2.top_bc * FT(83.26) +
                p.canopy.autotrophic_respiration.Ra # [3.664 kg CO2/ kg C] x [10^3 g CO2/ kg CO2] x [1 mol CO2/44.009 g CO2] = 83.26 mol CO2/kg C
        end
        return out
    else
        if isnothing(land_model.soilco2)
            @. out = p.canopy.autotrophic_respiration.Ra

        else
            out .=
                p.soilco2.top_bc .* FT(83.26) .+
                p.canopy.autotrophic_respiration.Ra
        end
    end
end
@diagnostic_compute "total_respiration" CanopyModel p.canopy.autotrophic_respiration.Ra

function compute_latent_heat_flux!(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel{FT},
) where {FT}
    if isnothing(out)
        out = zeros(land_model.canopy.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        @. out = p.soil.turbulent_fluxes.lhf + p.canopy.turbulent_fluxes.lhf
        return out
    else
        out .= p.soil.turbulent_fluxes.lhf .+ p.canopy.turbulent_fluxes.lhf
    end
end

function compute_sensible_heat_flux!(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel{FT},
) where {FT}
    if isnothing(out)
        out = zeros(land_model.canopy.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        @. out = p.soil.turbulent_fluxes.shf + p.canopy.turbulent_fluxes.shf
        return out
    else
        out .= p.soil.turbulent_fluxes.shf .+ p.canopy.turbulent_fluxes.shf
    end
end

function compute_latent_heat_flux!(
    out,
    Y,
    p,
    t,
    land_model::LandModel{FT},
) where {FT}
    if isnothing(out)
        out = zeros(land_model.canopy.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        @. out =
            p.soil.turbulent_fluxes.lhf * (1 - p.snow.snow_cover_fraction) +
            p.canopy.turbulent_fluxes.lhf +
            p.snow.snow_cover_fraction * p.snow.turbulent_fluxes.lhf
        return out
    else
        @. out =
            p.soil.turbulent_fluxes.lhf * (1 - p.snow.snow_cover_fraction) +
            p.canopy.turbulent_fluxes.lhf +
            p.snow.snow_cover_fraction * p.snow.turbulent_fluxes.lhf
    end
end

function compute_sensible_heat_flux!(
    out,
    Y,
    p,
    t,
    land_model::LandModel{FT},
) where {FT}
    if isnothing(out)
        out = zeros(land_model.canopy.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        @. out =
            p.soil.turbulent_fluxes.shf * (1 - p.snow.snow_cover_fraction) +
            p.canopy.turbulent_fluxes.shf +
            p.snow.snow_cover_fraction * p.snow.turbulent_fluxes.shf
        return out
    else
        @. out =
            p.soil.turbulent_fluxes.shf * (1 - p.snow.snow_cover_fraction) +
            p.canopy.turbulent_fluxes.shf +
            p.snow.snow_cover_fraction * p.snow.turbulent_fluxes.shf
    end
end

function compute_net_radiation!(
    out,
    Y,
    p,
    t,
    land_model::Union{SoilCanopyModel{FT}, LandModel{FT}},
) where {FT}
    if isnothing(out)
        out = zeros(land_model.canopy.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        @. out = p.drivers.LW_d - p.LW_u + p.drivers.SW_d - p.SW_u
        return out

    else
        out .= p.drivers.LW_d .- p.LW_u .+ p.drivers.SW_d .- p.SW_u

    end
end

# variables stored in Y (prognostic or state variables)
nan_if_no_canopy(T::FT, PAI::FT) where {FT <: Real} = PAI > 0 ? T : FT(NaN)
function compute_canopy_temperature!(
    out,
    Y,
    p,
    t,
    land_model::Union{SoilCanopyModel{FT}, LandModel{FT}},
) where {FT}
    PAI = p.scratch1
    @. PAI = p.canopy.biomass.area_index.leaf + p.canopy.biomass.area_index.stem
    if isnothing(out)
        out = zeros(land_model.canopy.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        out .=
            nan_if_no_canopy.(
                canopy_temperature(land_model.canopy.energy, land_model, Y, p),
                PAI,
            )
        return out
    else
        out .=
            nan_if_no_canopy.(
                canopy_temperature(land_model.canopy.energy, land_model, Y, p),
                PAI,
            )
    end
end
function compute_canopy_temperature!(
    out,
    Y,
    p,
    t,
    land_model::CanopyModel{FT},
) where {FT}
    PAI = p.canopy.biomass.area_index.leaf .+ p.canopy.biomass.area_index.stem # Allocates
    if isnothing(out)
        out = zeros(land_model.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        out .=
            nan_if_no_canopy.(
                canopy_temperature(land_model.energy, land_model, Y, p),
                PAI,
            )
        return out
    else
        out .=
            nan_if_no_canopy.(
                canopy_temperature(land_model.energy, land_model, Y, p),
                PAI,
            )
    end
end

@diagnostic_compute "soilco2" Union{SoilCanopyModel, LandModel} Y.soilco2.CO2
@diagnostic_compute "soilo2" Union{SoilCanopyModel, LandModel} Y.soilco2.O2_f
@diagnostic_compute "soc" Union{SoilCanopyModel, LandModel} Y.soilco2.SOC
@diagnostic_compute "soil_water_content" Union{
    SoilCanopyModel,
    LandModel,
    EnergyHydrology,
} Y.soil.ϑ_l
@diagnostic_compute "soil_ice_content" Union{
    SoilCanopyModel,
    LandModel,
    EnergyHydrology,
} Y.soil.θ_i
@diagnostic_compute "soil_internal_energy" Union{
    SoilCanopyModel,
    LandModel,
    EnergyHydrology,
} Y.soil.ρe_int

@diagnostic_compute "snow_water_equivalent" LandModel Y.snow.S
@diagnostic_compute "snow_depth" LandModel p.snow.z_snow
@diagnostic_compute "snow_cover_fraction" LandModel p.snow.snow_cover_fraction
@diagnostic_compute "evapotranspiration" EnergyHydrology p.soil.turbulent_fluxes.vapor_flux_liq

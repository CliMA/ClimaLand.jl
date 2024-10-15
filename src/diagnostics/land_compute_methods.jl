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
        return copy(p.soil.R_n)
    else
        out .= p.soil.R_n
    end
end
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
                    return copy($compute)
                else
                    out .= $compute
                end
            end
        end,
    )
end



### BucketModel ###

# variables stored in p (diagnostics variables stored in the cache)
@diagnostic_compute "aerodynamic_resistance" BucketModel p.bucket.turbulent_fluxes.r_ae
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

### SoilCanopyModel ###

# variables stored in p (diagnostics variables stored in the cache)

## Canopy Module ##

# Canopy - Solar Induced Fluorescence
@diagnostic_compute "solar_induced_fluorescence" SoilCanopyModel p.canopy.sif.SIF

# Canopy - Autotrophic respiration
@diagnostic_compute "autotrophic_respiration" SoilCanopyModel p.canopy.autotrophic_respiration.Ra

# Canopy - Conductance
@diagnostic_compute "stomatal_conductance" SoilCanopyModel p.canopy.conductance.gs
@diagnostic_compute "canopy_transpiration" SoilCanopyModel p.canopy.energy.turbulent_fluxes.transpiration

# Canopy - Energy
@diagnostic_compute "canopy_latent_heat_flux" SoilCanopyModel p.canopy.energy.turbulent_fluxes.lhf
@diagnostic_compute "canopy_sensible_heat_flux" SoilCanopyModel p.canopy.energy.turbulent_fluxes.shf

# Canopy - Hydraulics
#@diagnostic_compute "leaf_water_potential" SoilCanopyModel last(
#    p.canopy.hydraulics.ψ,
#)
# @diagnostic_compute "flux_per_ground_area" SoilCanopyModel p.canopy.hydraulics.fa # return a Tuple
@diagnostic_compute "root_flux_per_ground_area" SoilCanopyModel p.canopy.hydraulics.fa_roots
@diagnostic_compute "leaf_area_index" SoilCanopyModel p.canopy.hydraulics.area_index.leaf
@diagnostic_compute "moisture_stress_factor" SoilCanopyModel p.canopy.hydraulics.β
@diagnostic_compute "root_area_index" SoilCanopyModel p.canopy.hydraulics.area_index.root
@diagnostic_compute "stem_area_index" SoilCanopyModel p.canopy.hydraulics.area_index.stem

# Canopy - Photosynthesis
@diagnostic_compute "photosynthesis_net_canopy" SoilCanopyModel p.canopy.photosynthesis.GPP
@diagnostic_compute "photosynthesis_net_leaf" SoilCanopyModel p.canopy.photosynthesis.An
@diagnostic_compute "respiration_leaf" SoilCanopyModel p.canopy.photosynthesis.Rd
@diagnostic_compute "vcmax25" SoilCanopyModel p.canopy.photosynthesis.Vcmax25

# Canopy - Radiative Transfer
@diagnostic_compute "near_infrared_radiation_down" SoilCanopyModel p.canopy.radiative_transfer.nir_d
@diagnostic_compute "near_infrared_radiation_absorbed" SoilCanopyModel p.canopy.radiative_transfer.nir.abs
@diagnostic_compute "near_infrared_radiation_reflected" SoilCanopyModel p.canopy.radiative_transfer.nir.refl
@diagnostic_compute "near_infrared_radiation_transmitted" SoilCanopyModel p.canopy.radiative_transfer.nir.trans
@diagnostic_compute "photosynthetically_active_radiation_down" SoilCanopyModel p.canopy.radiative_transfer.par_d
@diagnostic_compute "photosynthetically_active_radiation_absorbed" SoilCanopyModel p.canopy.radiative_transfer.par.abs
@diagnostic_compute "photosynthetically_active_radiation_reflected" SoilCanopyModel p.canopy.radiative_transfer.par.refl
@diagnostic_compute "photosynthetically_active_radiation_transmitted" SoilCanopyModel p.canopy.radiative_transfer.par.trans
@diagnostic_compute "radiation_longwave_net" SoilCanopyModel p.canopy.radiative_transfer.LW_n
@diagnostic_compute "radiation_shortwave_net" SoilCanopyModel p.canopy.radiative_transfer.SW_n

## Drivers Module ##

@diagnostic_compute "soil_organic_carbon" SoilCanopyModel p.drivers.soc # need to fix this in src/shared_utilities/drivers
@diagnostic_compute "pressure" SoilCanopyModel p.drivers.P
@diagnostic_compute "rainfall" SoilCanopyModel p.drivers.P_liq
@diagnostic_compute "radiation_longwave_down" SoilCanopyModel p.drivers.LW_d
@diagnostic_compute "radiation_shortwave_down" SoilCanopyModel p.drivers.SW_d
@diagnostic_compute "snowfall" SoilCanopyModel p.drivers.P_snow
@diagnostic_compute "solar_zenith_angle" SoilCanopyModel p.drivers.θs
@diagnostic_compute "specific_humidity" SoilCanopyModel p.drivers.q
@diagnostic_compute "wind_speed" SoilCanopyModel p.drivers.u

## Soil Module ##

@diagnostic_compute "infiltration" SoilCanopyModel p.soil.infiltration
@diagnostic_compute "soil_hydraulic_conductivity" SoilCanopyModel p.soil.K
@diagnostic_compute "soil_thermal_conductivity" SoilCanopyModel p.soil.κ
@diagnostic_compute "soil_water_potential" SoilCanopyModel p.soil.ψ
@diagnostic_compute "soil_net_radiation" SoilCanopyModel p.soil.R_n
@diagnostic_compute "soil_temperature" SoilCanopyModel p.soil.T

# Soil - Turbulent Fluxes
@diagnostic_compute "soil_latent_heat_flux" SoilCanopyModel p.soil.turbulent_fluxes.lhf
@diagnostic_compute "soil_sensible_heat_flux" SoilCanopyModel p.soil.turbulent_fluxes.shf
@diagnostic_compute "vapor_flux" SoilCanopyModel p.soil.turbulent_fluxes.vapor_flux_liq # should add ice here

# Soil - SoilCO2
function compute_heterotrophic_respiration!(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel{FT},
) where {FT}
    if isnothing(out)
        return p.soilco2.top_bc .* FT(83.26)
    else
        out .= p.soilco2.top_bc .* FT(83.26)
    end
end # Convert from kg C to mol CO2.
# To convert from kg C to mol CO2, we need to multiply by:
# [3.664 kg CO2/ kg C] x [10^3 g CO2/ kg CO2] x [1 mol CO2/44.009 g CO2] = 83.26 mol CO2/kg C

@diagnostic_compute "soilco2_diffusivity" SoilCanopyModel p.soilco2.D
@diagnostic_compute "soilco2_source_microbe" SoilCanopyModel p.soilco2.Sm

## Other ##
@diagnostic_compute "sw_albedo" SoilCanopyModel p.α_sfc
@diagnostic_compute "lw_up" SoilCanopyModel p.LW_u
@diagnostic_compute "sw_up" SoilCanopyModel p.SW_u
@diagnostic_compute "surface_runoff" SoilCanopyModel p.soil.R_s

function compute_evapotranspiration!(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel{FT},
) where {FT}
    if isnothing(out)
        return (
            p.soil.turbulent_fluxes.vapor_flux_liq .+
            p.soil.turbulent_fluxes.vapor_flux_ice .+
            p.canopy.energy.turbulent_fluxes.transpiration
        ) .* 1000 # density of liquid water (1000kg/m^3)
    else
        out .=
            (
                p.soil.turbulent_fluxes.vapor_flux_liq .+
                p.soil.turbulent_fluxes.vapor_flux_ice .+
                p.canopy.energy.turbulent_fluxes.transpiration
            ) .* 1000 # density of liquid water (1000kg/m^3)
    end
end

function compute_total_respiration!(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel{FT},
) where {FT}
    if isnothing(out)
        return p.soilco2.top_bc .* FT(83.26) .+ # [3.664 kg CO2/ kg C] x [10^3 g CO2/ kg CO2] x [1 mol CO2/44.009 g CO2] = 83.26 mol CO2/kg C
               p.canopy.autotrophic_respiration.Ra
    else
        out .=
            p.soilco2.top_bc .* FT(83.26) .+ p.canopy.autotrophic_respiration.Ra
    end
end

# variables stored in Y (prognostic or state variables)
nan_if_no_canopy(T::FT, AI::FT) where {FT <: Real} = AI > 0 ? T : FT(NaN)
function compute_canopy_temperature!(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel{FT},
) where {FT}
    AI =
        p.canopy.hydraulics.area_index.leaf .+
        p.canopy.hydraulics.area_index.stem
    if isnothing(out)
        return nan_if_no_canopy.(Y.canopy.energy.T, AI)
    else
        out .= nan_if_no_canopy.(Y.canopy.energy.T, AI)
    end
end

@diagnostic_compute "soilco2" SoilCanopyModel Y.soilco2.C
@diagnostic_compute "soil_water_content" SoilCanopyModel Y.soil.ϑ_l
# @diagnostic_compute "plant_water_content" SoilCanopyModel Y.canopy.hydraulics.ϑ_l # return a Tuple
@diagnostic_compute "soil_ice_content" SoilCanopyModel Y.soil.θ_i
@diagnostic_compute "soil_internal_energy" SoilCanopyModel Y.soil.ρe_int

### EnergyHydrology ###

@diagnostic_compute "soil_water_content" EnergyHydrology Y.soil.ϑ_l
@diagnostic_compute "soil_ice_content" EnergyHydrology Y.soil.θ_i
@diagnostic_compute "soil_internal_energy" EnergyHydrology Y.soil.ρe_int
@diagnostic_compute "soil_temperature" EnergyHydrology p.soil.T

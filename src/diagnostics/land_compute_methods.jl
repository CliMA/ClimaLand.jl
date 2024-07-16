"""
    @land_compute name model compute


Macro generating a function to compute a land diagnostic,
needed in the ClimaDiagnostics framework.

For example,
@land_compute "soil_net_radiation" SoilCanopyModel p.soil.R_n

generates the function

function compute_soil_net_radiation!(out, Y, p, t, land_model::SoilCanopyModel)
    if isnothing(out)
        return copy(p.soil.R_n)
    else
        out .= p.soil.R_n
    end
end
"""
macro land_compute(name, model, compute)
    functions = Expr[]
    function_name = Symbol("compute_", name, "!")
    function_body = esc(quote
        function $function_name(out, Y, p, t, land_model::$model)
            if isnothing(out)
                return copy($compute)
            else
                out .= $compute
            end
        end
    end)
    push!(functions, function_body)
    return quote
        $(functions...)
    end
end

### BucketModel ###

# variables stored in p (diagnostics variables stored in the cache)
@land_compute "aerodynamic_resistance" BucketModel p.bucket.turbulent_fluxes.r_ae
@land_compute "albedo" BucketModel p.bucket.α_sfc
@land_compute "latent_heat_flux" BucketModel p.bucket.turbulent_fluxes.lhf
@land_compute "net_radiation" BucketModel p.bucket.R_n
@land_compute "sensible_heat_flux" BucketModel p.bucket.turbulent_fluxes.shf
@land_compute "surface_air_density" BucketModel p.bucket.ρ_sfc
@land_compute "surface_specific_humidity" BucketModel p.bucket.q_sfc
@land_compute "surface_temperature" BucketModel p.bucket.T_sfc
@land_compute "vapor_flux" BucketModel p.bucket.turbulent_fluxes.vapor_flux

# variables stored in Y (prognostic or state variables)
@land_compute "snow_water_equivalent" BucketModel Y.bucket.σS
@land_compute "soil_temperature" BucketModel Y.bucket.T
@land_compute "subsurface_water_storage" BucketModel Y.bucket.W
@land_compute "surface_water_content" BucketModel Y.bucket.Ws

### SoilCanopyModel ###

# variables stored in p (diagnostics variables stored in the cache)

## Canopy Module ##

# Canopy - Solar Induced Fluorescence
@land_compute "solar_induced_fluorescence" SoilCanopyModel p.canopy.sif.SIF

# Canopy - Autotrophic respiration
@land_compute "autotrophic_respiration" SoilCanopyModel p.canopy.autotrophic_respiration.Ra

# Canopy - Conductance
@land_compute "stomatal_conductance" SoilCanopyModel p.canopy.conductance.gs
@land_compute "canopy_transpiration" SoilCanopyModel p.canopy.conductance.transpiration

# Canopy - Energy
@land_compute "canopy_aerodynamic_resistance" SoilCanopyModel p.canopy.energy.r_ae
@land_compute "canopy_latent_heat_flux" SoilCanopyModel p.canopy.energy.lhf
@land_compute "canopy_sensible_heat_flux" SoilCanopyModel p.canopy.energy.shf

# Canopy - Hydraulics
@land_compute "canopy_water_potential" SoilCanopyModel p.canopy.hydraulics.ψ
@land_compute "cross_section" SoilCanopyModel p.canopy.hydraulics.fa
@land_compute "cross_section_roots" SoilCanopyModel p.canopy.hydraulics.fa_roots
@land_compute "leaf_area_index" SoilCanopyModel p.canopy.hydraulics.area_index.leaf
@land_compute "moisture_stress_factor" SoilCanopyModel p.canopy.hydraulics.β
@land_compute "root_area_index" SoilCanopyModel p.canopy.hydraulics.area_index.root
@land_compute "stem_area_index" SoilCanopyModel p.canopy.hydraulics.area_index.stem

# Canopy - Photosynthesis
@land_compute "photosynthesis_net_canopy" SoilCanopyModel p.canopy.photosynthesis.GPP
@land_compute "photosynthesis_net_leaf" SoilCanopyModel p.canopy.photosynthesis.An
@land_compute "respiration_leaf" SoilCanopyModel p.canopy.photosynthesis.Rd
@land_compute "vcmax25" SoilCanopyModel p.canopy.photosynthesis.Vcmax25

# Canopy - Radiative Transfer
@land_compute "near_infrared_radiation" SoilCanopyModel p.canopy.radiative_transfer.inc_nir
@land_compute "near_infrared_radiation_absorbed" SoilCanopyModel p.canopy.radiative_transfer.nir.abs
@land_compute "near_infrared_radiation_reflected" SoilCanopyModel p.canopy.radiative_transfer.nir.refl
@land_compute "near_infrared_radiation_transmitted" SoilCanopyModel p.canopy.radiative_transfer.nir.trans
@land_compute "photosynthetically_active_radiation" SoilCanopyModel p.canopy.radiative_transfer.inc_par
@land_compute "photosynthetically_active_radiation_absorbed" SoilCanopyModel p.canopy.radiative_transfer.par.abs
@land_compute "photosynthetically_active_radiation_reflected" SoilCanopyModel p.canopy.radiative_transfer.par.refl
@land_compute "photosynthetically_active_radiation_transmitted" SoilCanopyModel p.canopy.radiative_transfer.par.trans
@land_compute "radiation_longwave_net" SoilCanopyModel p.canopy.radiative_transfer.LW_n
@land_compute "radiation_shortwave_net" SoilCanopyModel p.canopy.radiative_transfer.SW_n

## Drivers Module ##

@land_compute "soil_organic_carbon" SoilCanopyModel p.drivers.soc # need to fix this in src/shared_utilities/drivers 
@land_compute "pressure" SoilCanopyModel p.drivers.P
@land_compute "rainfall" SoilCanopyModel p.drivers.P_liq
@land_compute "radiation_longwave_down" SoilCanopyModel p.drivers.LW_d
@land_compute "radiation_shortwave_down" SoilCanopyModel p.drivers.SW_d
@land_compute "snowfall" SoilCanopyModel p.drivers.P_snow
@land_compute "solar_zenith_angle" SoilCanopyModel p.drivers.θs
@land_compute "specific_humidity" SoilCanopyModel p.drivers.q
@land_compute "wind_speed" SoilCanopyModel p.drivers.u

## Soil Module ##

@land_compute "infiltration" SoilCanopyModel p.soil.infiltration
@land_compute "soil_hydraulic_conductivity" SoilCanopyModel p.soil.K
@land_compute "soil_thermal_conductivity" SoilCanopyModel p.soil.κ
@land_compute "soil_water_potential" SoilCanopyModel p.soil.ψ
@land_compute "soil_net_radiation" SoilCanopyModel p.soil.R_n
@land_compute "soil_temperature" SoilCanopyModel p.soil.T

# Soil - Turbulent Fluxes
@land_compute "soil_aerodynamic_resistance" SoilCanopyModel p.soil.turbulent_fluxes.r_ae
@land_compute "soil_latent_heat_flux" SoilCanopyModel p.soil.turbulent_flux.lhf
@land_compute "soil_sensible_heat_flux" SoilCanopyModel p.soil.turbulent_fluxes.shf
@land_compute "vapor_flux" SoilCanopyModel p.soil.turbulent_fluxes.vapor_flux

# Soil - SoilCO2
@land_compute "heterotrophic_respiration" SoilCanopyModel p.soilco2.top_bc
@land_compute "soilco2_diffusivity" SoilCanopyModel p.soilco2.D
@land_compute "soilco2_source_microbe" SoilCanopyModel p.soilco2.Sm

# variables stored in Y (prognostic or state variables)

@land_compute "canopy_temperature" SoilCanopyModel Y.canopy.energy.T
@land_compute "soilco2" SoilCanopyModel Y.soilco2.C
@land_compute "soil_water_content" SoilCanopyModel Y.soil.ϑ_l
@land_compute "plant_water_content" SoilCanopyModel Y.canopy.hydraulics.ϑ_l
@land_compute "soil_ice_content" SoilCanopyModel Y.soil.θ_i
@land_compute "soil_internal_energy" SoilCanopyModel Y.soil.ρe_int

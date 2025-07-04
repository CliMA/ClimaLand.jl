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

### Conservation ##
@diagnostic_compute "water_volume_per_area" EnergyHydrology p.soil.total_water
@diagnostic_compute "energy_per_area" EnergyHydrology p.soil.total_energy
@diagnostic_compute "water_volume_per_area_change" EnergyHydrology Y.soil.∫F_vol_liq_water_dt
@diagnostic_compute "energy_per_area_change" EnergyHydrology Y.soil.∫F_e_dt


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

### Union{SoilCanopyModel, LandModel} ###

# variables stored in p (diagnostics variables stored in the cache)

## Canopy Module ##

# Canopy - Solar Induced Fluorescence
@diagnostic_compute "solar_induced_fluorescence" Union{
    SoilCanopyModel,
    LandModel,
} p.canopy.sif.SIF

# Canopy - Autotrophic respiration
@diagnostic_compute "autotrophic_respiration" Union{SoilCanopyModel, LandModel} p.canopy.autotrophic_respiration.Ra

# Canopy - Conductance
function compute_stomatal_conductance!(
    out,
    Y,
    p,
    t,
    land_model::Union{SoilCanopyModel, LandModel},
)
    canopy = land_model.canopy
    (; g1, g0, Drel) = canopy.conductance.parameters
    earth_param_set = canopy.parameters.earth_param_set
    thermo_params = LP.thermodynamic_parameters(earth_param_set)

    if isnothing(out)
        out = zeros(land_model.canopy.domain.space.surface) # Allocates
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
            p.canopy.photosynthesis.An,
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
            p.canopy.photosynthesis.An,
            p.drivers.c_co2,
        )
    end
end

function compute_canopy_transpiration!(
    out,
    Y,
    p,
    t,
    land_model::Union{SoilCanopyModel{FT}, LandModel{FT}},
) where {FT}
    # Convert to a mass flux by multiplying by the density of liquid
    # water
    if isnothing(out)
        out = zeros(land_model.canopy.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        @. out = p.canopy.turbulent_fluxes.transpiration * 1000
        return out
    else
        @. out = p.canopy.turbulent_fluxes.transpiration * 1000
    end
end

# Canopy - Energy
@diagnostic_compute "canopy_latent_heat_flux" Union{SoilCanopyModel, LandModel} p.canopy.turbulent_fluxes.lhf
@diagnostic_compute "canopy_sensible_heat_flux" Union{
    SoilCanopyModel,
    LandModel,
} p.canopy.turbulent_fluxes.shf

# Canopy - Hydraulics
function compute_leaf_water_potential!(
    out,
    Y,
    p,
    t,
    land_model::Union{SoilCanopyModel, LandModel},
)
    hydraulics = land_model.canopy.hydraulics
    n_stem = hydraulics.n_stem
    n_leaf = hydraulics.n_leaf
    n = n_stem + n_leaf
    if isnothing(out)
        out = zeros(land_model.canopy.domain.space.surface) # Allocates
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
} p.canopy.hydraulics.fa_roots
@diagnostic_compute "leaf_area_index" Union{SoilCanopyModel, LandModel} p.canopy.hydraulics.area_index.leaf
# Canopy - Hydraulics
function compute_moisture_stress_factor!(
    out,
    Y,
    p,
    t,
    land_model::Union{SoilCanopyModel, LandModel},
)
    canopy = land_model.canopy

    hydraulics = canopy.hydraulics
    n_stem = hydraulics.n_stem
    n_leaf = hydraulics.n_leaf
    n = n_stem + n_leaf

    earth_param_set = canopy.parameters.earth_param_set
    grav = LP.grav(earth_param_set)
    ρ_l = LP.ρ_cloud_liq(earth_param_set)
    (; sc, pc) = canopy.photosynthesis.parameters
    ψ = p.canopy.hydraulics.ψ
    if isnothing(out)
        out = zeros(land_model.canopy.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        @. out = moisture_stress(ψ.:($$n) * ρ_l * grav, sc, pc)
        return out
    else
        @. out = moisture_stress(ψ.:($$n) * ρ_l * grav, sc, pc)
    end
end
@diagnostic_compute "root_area_index" Union{SoilCanopyModel, LandModel} p.canopy.hydraulics.area_index.root
@diagnostic_compute "stem_area_index" Union{SoilCanopyModel, LandModel} p.canopy.hydraulics.area_index.stem

# Canopy - Photosynthesis
@diagnostic_compute "photosynthesis_net_canopy" Union{
    SoilCanopyModel,
    LandModel,
} p.canopy.photosynthesis.GPP
@diagnostic_compute "photosynthesis_net_leaf" Union{SoilCanopyModel, LandModel} p.canopy.photosynthesis.An
@diagnostic_compute "respiration_leaf" Union{SoilCanopyModel, LandModel} p.canopy.photosynthesis.Rd

# Canopy - Radiative Transfer
@diagnostic_compute "near_infrared_radiation_down" Union{
    SoilCanopyModel,
    LandModel,
} p.canopy.radiative_transfer.nir_d
@diagnostic_compute "near_infrared_radiation_absorbed" Union{
    SoilCanopyModel,
    LandModel,
} p.canopy.radiative_transfer.nir.abs
@diagnostic_compute "near_infrared_radiation_reflected" Union{
    SoilCanopyModel,
    LandModel,
} p.canopy.radiative_transfer.nir.refl
@diagnostic_compute "near_infrared_radiation_transmitted" Union{
    SoilCanopyModel,
    LandModel,
} p.canopy.radiative_transfer.nir.trans
@diagnostic_compute "photosynthetically_active_radiation_down" Union{
    SoilCanopyModel,
    LandModel,
} p.canopy.radiative_transfer.par_d
@diagnostic_compute "photosynthetically_active_radiation_absorbed" Union{
    SoilCanopyModel,
    LandModel,
} p.canopy.radiative_transfer.par.abs
@diagnostic_compute "photosynthetically_active_radiation_reflected" Union{
    SoilCanopyModel,
    LandModel,
} p.canopy.radiative_transfer.par.refl
@diagnostic_compute "photosynthetically_active_radiation_transmitted" Union{
    SoilCanopyModel,
    LandModel,
} p.canopy.radiative_transfer.par.trans
@diagnostic_compute "radiation_longwave_net" Union{SoilCanopyModel, LandModel} p.canopy.radiative_transfer.LW_n
@diagnostic_compute "radiation_shortwave_net" Union{SoilCanopyModel, LandModel} p.canopy.radiative_transfer.SW_n

## Drivers Module ##

@diagnostic_compute "soil_organic_carbon" Union{SoilCanopyModel, LandModel} p.drivers.soc # need to fix this in src/shared_utilities/drivers
@diagnostic_compute "pressure" Union{SoilCanopyModel, LandModel} p.drivers.P
@diagnostic_compute "rainfall" Union{SoilCanopyModel, LandModel} p.drivers.P_liq
@diagnostic_compute "radiation_longwave_down" Union{SoilCanopyModel, LandModel} p.drivers.LW_d
@diagnostic_compute "radiation_shortwave_down" Union{SoilCanopyModel, LandModel} p.drivers.SW_d
@diagnostic_compute "snowfall" Union{SoilCanopyModel, LandModel} p.drivers.P_snow
@diagnostic_compute "tair" Union{SoilCanopyModel, LandModel} p.drivers.T
@diagnostic_compute "specific_humidity" Union{SoilCanopyModel, LandModel} p.drivers.q
@diagnostic_compute "wind_speed" Union{SoilCanopyModel, LandModel} p.drivers.u

function compute_precip!(
    out,
    Y,
    p,
    t,
    land_model::Union{SoilCanopyModel{FT}, LandModel{FT}},
) where {FT}
    if isnothing(out)
        out = zeros(land_model.soil.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        @. out = (p.drivers.P_liq + p.drivers.P_snow) * 1000 # density of liquid water (1000kg/m^3)
        return out
    else
        @. out = (p.drivers.P_liq + p.drivers.P_snow) * 1000# density of liquid water (1000kg/m^3)
    end
end

## Soil Module ##

@diagnostic_compute "infiltration" Union{SoilCanopyModel, LandModel} p.soil.infiltration
@diagnostic_compute "soil_hydraulic_conductivity" Union{
    SoilCanopyModel,
    LandModel,
} p.soil.K
@diagnostic_compute "soil_thermal_conductivity" Union{
    SoilCanopyModel,
    LandModel,
} p.soil.κ
@diagnostic_compute "soil_water_potential" Union{SoilCanopyModel, LandModel} p.soil.ψ
@diagnostic_compute "soil_net_radiation" Union{SoilCanopyModel, LandModel} p.soil.R_n
@diagnostic_compute "soil_temperature" Union{SoilCanopyModel, LandModel} p.soil.T

function compute_10cm_water_mass!(
    out,
    Y,
    p,
    t,
    land_model::Union{SoilCanopyModel{FT}, LandModel{FT}},
) where {FT}
    ∫Hθdz = p.soil.sfc_scratch
    Hθ = p.soil.sub_sfc_scratch
    z = land_model.soil.domain.fields.z
    depth = FT(-0.1)
    earth_param_set = land_model.soil.parameters.earth_param_set
    _ρ_liq = LP.ρ_cloud_liq(earth_param_set)
    _ρ_ice = LP.ρ_cloud_ice(earth_param_set)
    # Convert from volumetric water content to water mass per unit volume using density
    @. Hθ = (Y.soil.ϑ_l * _ρ_liq + Y.soil.θ_i * _ρ_ice) * heaviside(z, depth)
    column_integral_definite!(∫Hθdz, Hθ)

    # The layering of the soil model may not coincide with 10 cm exactly, and this could lead
    # to the integral above not exactly representing 10cm.
    # To adjust, divide by the ∫heaviside(z, depth) dz, and then multiply by 10cm
    H = p.subsfc_scratch
    @. H = heaviside(z, depth)
    ∫Hdz = p.sfc_scratch
    column_integral_definite!(∫Hdz, H)

    if isnothing(out)
        out = zeros(land_model.soil.domain.space.surface) # Allocates
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
    land_model::SoilCanopyModel{FT},
) where {FT}
    if isnothing(out)
        out = zeros(land_model.soil.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        @. out = (p.soil.PAR_albedo + p.soil.NIR_albedo) / 2
        return out
    else
        @. out = (p.soil.PAR_albedo + p.soil.NIR_albedo) / 2
    end
end

# Soil - Turbulent Fluxes
@diagnostic_compute "soil_latent_heat_flux" Union{SoilCanopyModel, LandModel} p.soil.turbulent_fluxes.lhf
@diagnostic_compute "soil_sensible_heat_flux" Union{SoilCanopyModel, LandModel} p.soil.turbulent_fluxes.shf
@diagnostic_compute "vapor_flux" Union{SoilCanopyModel, LandModel} p.soil.turbulent_fluxes.vapor_flux_liq # should add ice here

# Soil - SoilCO2
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

@diagnostic_compute "soilco2_diffusivity" Union{SoilCanopyModel, LandModel} p.soilco2.D
@diagnostic_compute "soilco2_source_microbe" Union{SoilCanopyModel, LandModel} p.soilco2.Sm

## Other ##
@diagnostic_compute "sw_albedo" Union{SoilCanopyModel, LandModel} p.α_sfc
@diagnostic_compute "lw_up" Union{SoilCanopyModel, LandModel} p.LW_u
@diagnostic_compute "sw_up" Union{SoilCanopyModel, LandModel} p.SW_u
@diagnostic_compute "surface_runoff" Union{SoilCanopyModel, LandModel} p.soil.R_s
@diagnostic_compute "subsurface_runoff" Union{SoilCanopyModel, LandModel} p.soil.R_ss

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
                p.canopy.turbulent_fluxes.transpiration
            ) * 1000 # density of liquid water (1000kg/m^3)
        return out
    else
        out .=
            (
                p.soil.turbulent_fluxes.vapor_flux_liq .+
                p.soil.turbulent_fluxes.vapor_flux_ice .+
                p.canopy.turbulent_fluxes.transpiration
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
                p.canopy.turbulent_fluxes.transpiration +
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
                p.canopy.turbulent_fluxes.transpiration +
                p.snow.snow_cover_fraction * p.snow.turbulent_fluxes.vapor_flux
            ) * 1000 # density of liquid water (1000kg/m^3)
    end
end

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
        @. out =
            p.soilco2.top_bc * FT(83.26) + p.canopy.autotrophic_respiration.Ra # [3.664 kg CO2/ kg C] x [10^3 g CO2/ kg CO2] x [1 mol CO2/44.009 g CO2] = 83.26 mol CO2/kg C
        return out
    else
        out .=
            p.soilco2.top_bc .* FT(83.26) .+ p.canopy.autotrophic_respiration.Ra
    end
end

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
nan_if_no_canopy(T::FT, AI::FT) where {FT <: Real} = AI > 0 ? T : FT(NaN)
function compute_canopy_temperature!(
    out,
    Y,
    p,
    t,
    land_model::Union{SoilCanopyModel{FT}, LandModel{FT}},
) where {FT}
    AI = p.scratch1
    @. AI =
        p.canopy.hydraulics.area_index.leaf +
        p.canopy.hydraulics.area_index.stem
    if isnothing(out)
        out = zeros(land_model.canopy.domain.space.surface) # Allocates
        fill!(field_values(out), NaN) # fill with NaNs, even over the ocean
        @. out = nan_if_no_canopy(Y.canopy.energy.T, AI)
        return out
    else
        out .= nan_if_no_canopy.(Y.canopy.energy.T, AI)
    end
end

@diagnostic_compute "soilco2" Union{SoilCanopyModel, LandModel} Y.soilco2.C
@diagnostic_compute "soil_water_content" Union{SoilCanopyModel, LandModel} Y.soil.ϑ_l
# @diagnostic_compute "plant_water_content" Union{SoilCanopyModel, LandModel} Y.canopy.hydraulics.ϑ_l # return a Tuple
@diagnostic_compute "soil_ice_content" Union{SoilCanopyModel, LandModel} Y.soil.θ_i
@diagnostic_compute "soil_internal_energy" Union{SoilCanopyModel, LandModel} Y.soil.ρe_int

@diagnostic_compute "snow_water_equivalent" LandModel Y.snow.S
@diagnostic_compute "snow_depth" LandModel p.snow.z_snow
@diagnostic_compute "snow_cover_fraction" LandModel p.snow.snow_cover_fraction

### EnergyHydrology ###

@diagnostic_compute "soil_water_content" EnergyHydrology Y.soil.ϑ_l
@diagnostic_compute "soil_ice_content" EnergyHydrology Y.soil.θ_i
@diagnostic_compute "soil_internal_energy" EnergyHydrology Y.soil.ρe_int
@diagnostic_compute "soil_temperature" EnergyHydrology p.soil.T
@diagnostic_compute "evapotranspiration" EnergyHydrology p.soil.turbulent_fluxes.vapor_flux_liq

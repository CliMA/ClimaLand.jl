# stored in p

function compute_soil_net_radiation!(out, Y, p, t, land_model::SoilCanopyModel)
    if isnothing(out)
        return copy(p.soil.R_n)
    else
        out .= p.soil.R_n
    end
end

function compute_soil_latent_heat_flux!(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.soil.turbulent_flux.lhf) # is this different from canopy.energy.lhf?
    else
        out .= p.soil.turbulent_flux.lhf
    end
end

function compute_soil_aerodynamic_resistance!(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.soil.turbulent_fluxes.r_ae)
    else
        out .= p.soil.turbulent_fluxes.r_ae
    end
end

function compute_soil_sensible_heat_flux!(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.soil.turbulent_fluxes.shf)
    else
        out .= p.soil.turbulent_fluxes.shf
    end
end

function compute_vapor_flux!(out, Y, p, t, land_model::SoilCanopyModel)
    if isnothing(out)
        return copy(p.soil.turbulent_fluxes.vapor_flux)
    else
        out .= p.soil.turbulent_fluxes.vapor_flux
    end
end

function compute_soil_temperature!(
    out,
    Y,
    p,
    t,
    land_model::Union{SoilCanopyModel, EnergyHydrology},
)
    if isnothing(out)
        return copy(top_center_to_surface(p.soil.T))
    else
        out .= top_center_to_surface(p.soil.T)
    end
end

function compute_soil_water_liquid!(
    out,
    Y,
    p,
    t,
    land_model::Union{SoilCanopyModel, EnergyHydrology},
)
    soil_params =
        typeof(land_model) <: SoilCanopyModel ? land_model.soil.parameters :
        land_model.parameters
    if isnothing(out)
        return copy(
            top_center_to_surface(
                (Y.soil.ϑ_l .- soil_params.θ_r) ./
                (soil_params.ν .- soil_params.θ_r),
            ),
        )
    else
        out .= top_center_to_surface(
            (Y.soil.ϑ_l .- soil_params.θ_r) ./
            (soil_params.ν .- soil_params.θ_r),
        )
    end
end

function compute_infiltration(
    out,
    Y,
    p,
    t,
    land_model::Union{SoilCanopyModel, EnergyHydrology},
)
    if isnothing(out)
        return copy(p.soil.infiltration)
    else
        out .= p.soil.infiltration
    end
end

function compute_soilco2_diffusivity(out, Y, p, t, land_model::SoilCanopyModel)
    if isnothing(out)
        return copy(p.soilco2.D) # NOTE: we will need a method to compute surface co2 efflux
    else
        out .= p.soilco2.D
    end
end

function compute_soilco2_source_microbe(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.soilco2.Sm)
    else
        out .= p.soilco2.Sm
    end
end

function compute_stomatal_conductance(out, Y, p, t, land_model::SoilCanopyModel)
    if isnothing(out)
        return copy(p.canopy.conductance.gs) # doublecheck: stomata, not canopy
    else
        out .= p.canopy.conductance.gs
    end
end

function compute_medlyn_term(out, Y, p, t, land_model::SoilCanopyModel)
    if isnothing(out)
        return copy(p.canopy.conductance.medlyn_term)
    else
        out .= p.canopy.conductance.medlyn_term
    end
end

function compute_canopy_transpiration(out, Y, p, t, land_model::SoilCanopyModel)
    if isnothing(out)
        return copy(p.canopy.transpiration) # doublecheck: canopy, not leaf
    else
        out .= p.canopy.transpiration
    end
end

function compute_rainfall(out, Y, p, t, land_model::SoilCanopyModel)
    if isnothing(out)
        return copy(p.drivers.P_liq) # I guess this is read and put it p. not computed. curious if we should handle this differently.
    else
        out .= p.drivers.P_liq
    end
end

function compute_snowfall(out, Y, p, t, land_model::SoilCanopyModel)
    if isnothing(out)
        return copy(p.drivers.P_snow) # following comment above, we could have a default getting only model output, and one also getting some inputs like drivers
    else
        out .= p.drivers.P_snow
    end
end

function compute_pressure(out, Y, p, t, land_model::SoilCanopyModel)
    if isnothing(out)
        return copy(p.drivers.P) # not sure if precip or pressure
    else
        out .= p.drivers.P
    end
end

function compute_wind_speed(out, Y, p, t, land_model::SoilCanopyModel)
    if isnothing(out)
        return copy(p.drivers.u)
    else
        out .= p.drivers.u
    end
end

function compute_specific_humidity(out, Y, p, t, land_model::SoilCanopyModel)
    if isnothing(out)
        return copy(p.drivers.q) # check if this is correct. Also, check if Bucket has it and same name or not.
    else
        out .= p.drivers.q
    end
end

function compute_air_co2(out, Y, p, t, land_model::SoilCanopyModel)
    if isnothing(out)
        return copy(p.drivers.c_co2)
    else
        out .= p.drivers.c_co2
    end
end

function compute_radiation_shortwave_down(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.drivers.SW_d)
    else
        out .= p.drivers.SW_d
    end
end

function compute_radiation_longwave_down(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.drivers.LW_d)
    else
        out .= p.drivers.LW_d
    end
end

function compute_photosynthesis_net_leaf(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.canopy.photosynthesis.An)
    else
        out .= p.canopy.photosynthesis.An
    end
end

function compute_photosynthesis_net_canopy!(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
) # could be gross primary productivity, but this is consistent with leaf
    if isnothing(out)
        return copy(p.canopy.photosynthesis.GPP)
    else
        out .= p.canopy.photosynthesis.GPP
    end
end

function compute_respiration_leaf(out, Y, p, t, land_model::SoilCanopyModel)
    if isnothing(out)
        return copy(p.canopy.photosynthesis.Rd)
    else
        out .= p.canopy.photosynthesis.Rd
    end
end

function compute_vcmax25(out, Y, p, t, land_model::SoilCanopyModel)
    if isnothing(out)
        return copy(p.canopy.photosynthesis.vcmax25)
    else
        out .= p.canopy.photosynthesis.vcmax25
    end
end

function compute_photosynthetically_active_radiation(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.canopy.radiative_transfer.par)
    else
        out .= p.canopy.radiative_transfer.par
    end
end

function compute_photosynthetically_active_radiation_absorbed(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.canopy.radiative_transfer.apar)
    else
        out .= p.canopy.radive_transfer.apar
    end
end

function compute_photosynthetically_active_radiation_reflected(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.canopy.radiative_transfer.rpar)
    else
        out .= p.canopy.radiative_transfer.rpar
    end
end

function compute_photosynthetically_active_radiation_transmitted(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.canopy.radiative_transfer.tpar)
    else
        out .= p.canopy.radiative_transfer.tpar
    end
end

function compute_near_infrared_radiation(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.canopy.radiative_transfer.nir)
    else
        out .= p.canopy.radiative_transfer.nir
    end
end

function compute_near_infrared_radiation_absorbed(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.canopy.radiative_transfer.anir)
    else
        out .= p.canopy.radiative_transfer.anir
    end
end

function compute_near_infrared_radiation_reflected(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.canopy.radiative_transfer.rnir)
    else
        out .= p.canopy.radiative_transfer.rnir
    end
end

function compute_near_infrared_radiation_transmitted(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.canopy.radiative_transfer.tnir)
    else
        out .= p.canopy.radiative_transfer.tnir
    end
end

function compute_radiation_shortwave_net(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.canopy.radiative_transfer.SW_n)
    else
        out .= p.canopy.radiative_transfer.SW_n
    end
end

function compute_radiation_longwave_net(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.canopy.radiative_transfer.LW_n)
    else
        out .= p.canopy.radiative_transfer.LW_n
    end
end

function compute_autotrophic_respiration(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.canopy.autotrophic_respiration.Ra)
    else
        out .= p.canopy.autotrophic_respiration.Ra
    end
end

function compute_soilco2(out, Y, p, t, land_model::SoilCanopyModel)
    if isnothing(out)
        return copy(Y.soilco2.C)
    else
        out .= Y.soilco2.C
    end
end

function compute_heterotrophic_respiration(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.soilco2.top_bc)
    else
        p.soilco2.top_bc
    end
end

function compute_soil_hydraulic_conductivity(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.soil.K)
    else
        p.soil.K
    end
end

function compute_soil_water_potential(out, Y, p, t, land_model::SoilCanopyModel)
    if isnothing(out)
        return copy(p.soil.ψ)
    else
        p.soil.ψ
    end
end

function compute_soil_thermal_conductivity(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.soil.κ)
    else
        p.soil.κ
    end
end

function compute_solar_zenith_angle(out, Y, p, t, land_model::SoilCanopyModel)
    if isnothing(out)
        return copy(p.drivers.θs)
    else
        p.drivers.θs
    end
end

function compute_moisture_stress_factor(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.canopy.hydraulics.β)
    else
        p.canopy.hydraulics.β
    end
end

function compute_canopy_water_potential(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.canopy.hydraulics.ψ)
    else
        p.canopy.hydraulics.ψ
    end
end

function compute_cross_section(out, Y, p, t, land_model::SoilCanopyModel)
    if isnothing(out)
        return copy(p.canopy.hydraulics.fa)
    else
        p.canopy.hydraulics.fa
    end
end

function compute_cross_section_roots(out, Y, p, t, land_model::SoilCanopyModel)
    if isnothing(out)
        return copy(p.canopy.hydraulics.fa_roots)
    else
        p.canopy.hydraulics.fa_roots
    end
end

function compute_area_index!(out, Y, p, t, land_model::SoilCanopyModel)
    if isnothing(out)
        return copy(p.canopy.hydraulics.area_index.leaf)
    else
        out .= p.canopy.hydraulics.area_index.leaf
    end
end

function compute_canopy_latent_heat_flux(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.canopy.energy.lhf)
    else
        p.canopy.canopy.lhf
    end
end

function compute_canopy_sensible_heat_flux(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.canopy.energy.shf)
    else
        p.canopy.canopy.shf
    end
end

function compute_canopy_aerodynamic_resistance(
    out,
    Y,
    p,
    t,
    land_model::SoilCanopyModel,
)
    if isnothing(out)
        return copy(p.canopy.energy.r_ae)
    else
        p.canopy.canopy.r_ae
    end
end

function compute_canopy_temperature!(out, Y, p, t, land_model::SoilCanopyModel)
    if isnothing(out)
        return copy(Y.canopy.energy.T)
    else
        out .= Y.canopy.energy.T
    end
end

function compute_soil_ice!(
    out,
    Y,
    p,
    t,
    land_model::Union{SoilCanopyModel, EnergyHydrology},
)
    soil_params =
        typeof(land_model) <: SoilCanopyModel ? land_model.soil.parameters :
        land_model.parameters
    if isnothing(out)
        return copy(
            top_center_to_surface(
                (Y.soil.θ_i .- soil_params.θ_r) ./
                (soil_params.ν .- soil_params.θ_r),
            ),
        )
    else
        out .= top_center_to_surface(
            (Y.soil.θ_i .- soil_params.θ_r) ./
            (soil_params.ν .- soil_params.θ_r),
        )
    end
end

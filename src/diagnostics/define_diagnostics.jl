# General helper functions for undefined diagnostics for a particular model
error_diagnostic_variable(variable, land_model::T) where {T} =
    error("Cannot compute $variable with model = $T")

# generate_error_functions is helper macro that generates the error message
# when the user tries calling something that is incompatible with the model
macro generate_error_functions(variable_names...)
    functions = Expr[]
    for variable in variable_names
        function_name_sym = Symbol("compute_", variable, "!")
        body = esc(quote
            function $function_name_sym(_, _, _, _, land_model)
                error_diagnostic_variable($variable, land_model)
            end
        end)
        push!(functions, body)
    end
    return quote
        $(functions...)
    end
end

# TODO: Automatically generate this list from the names of the diagnostics
@generate_error_functions "soil_net_radiation" "soil_latent_heat_flux" "soil_aerodynamic_resistance" "soil_sensible_heat_flux" "vapor_flux" "soil_temperature" "soil_water_liquid" "infiltration" "soilco2_diffusivity" "soilco2_source_microbe" "stomatal_conductance" "medlyn_term" "canopy_transpiration" "rainfall" "snowfall" "pressure" "wind_speed" "specific_humidity" "air_co2" "radiation_shortwave_down" "radiation_longwave_down" "photosynthesis_net_leaf" "photosynthesis_net_canopy" "respiration_leaf" "vcmax25" "photosynthetically_active_radiation" "photosynthetically_active_radiation_absorbed" "photosynthetically_active_radiation_reflected" "photosynthetically_active_radiation_transmitted" "near_infrared_radiation" "near_infrared_radiation_absorbed" "near_infrared_radiation_reflected" "near_infrared_radiation_transmitted" "radiation_shortwave_net" "radiation_longwave_net" "autotrophic_respiration" "soilco2" "heterotrophic_respiration" "soil_hydraulic_conductivity" "soil_water_potential" "soil_thermal_conductivity" "solar_zenith_angle" "moisture_stress_factor" "canopy_water_potential" "cross_section" "cross_section_roots" "area_index" "canopy_latent_heat_flux" "canopy_sensible_heat_flux" "canopy_aerodynamic_resistance" "canopy_temperature" "soil_ice"

"""
    define_diagnostics!(land_model)

Calls `add_diagnostic_variable!` for all available variables specializing the
compute function for `land_model`.
"""
function define_diagnostics!(land_model)

    # Stored in p

    # Albedo
    add_diagnostic_variable!(
        short_name = "alpha",
        long_name = "Albedo",
        standard_name = "albedo",
        units = "",
        comments = "The fraction of incoming radiation reflected by the land surface.",
        compute! = (out, Y, p, t) -> compute_albedo!(out, Y, p, t, land_model),
    )

    # Net radiation
    add_diagnostic_variable!(
        short_name = "rn",
        long_name = "Net Radiation",
        standard_name = "net_radiation",
        units = "W m^-2",
        comments = "Difference between incoming and outgoing shortwave and longwave radiation at the land surface.",
        compute! = (out, Y, p, t) ->
            compute_net_radiation!(out, Y, p, t, land_model),
    )

    # Surface temperature
    add_diagnostic_variable!(
        short_name = "tsfc",
        long_name = "Surface Temperature",
        standard_name = "surface_temperature",
        units = "K",
        comments = "Temperature of the land surface.",
        compute! = (out, Y, p, t) ->
            compute_surface_temperature!(out, Y, p, t, land_model),
    )

    # Surface specific humidity
    add_diagnostic_variable!(
        short_name = "qsfc",
        long_name = "Surface Specific Humidity",
        standard_name = "surface_specific_humidity",
        units = "",
        comments = "Ratio of water vapor mass to total moist air parcel mass.",
        compute! = (out, Y, p, t) ->
            compute_surface_specific_humidity!(out, Y, p, t, land_model),
    )

    # Latent heat flux
    add_diagnostic_variable!(
        short_name = "lhf",
        long_name = "Latent Heat Flux",
        standard_name = "latent_heat_flux",
        units = "W m^-2",
        comments = "Exchange of energy at the land-atmosphere interface due to water evaporation or sublimation.",
        compute! = (out, Y, p, t) ->
            compute_latent_heat_flux!(out, Y, p, t, land_model),
    )

    # Aerodynamic resistance
    add_diagnostic_variable!(
        short_name = "rae",
        long_name = "Aerodynamic Resistance",
        standard_name = "aerodynamic_resistance",
        units = "m s^-1",
        comments = "Effiency of turbulent transport controlling the land-atmosphere exchange of sensible and latent heat.",
        compute! = (out, Y, p, t) ->
            compute_aerodynamic_resistance!(out, Y, p, t, land_model),
    )

    # Sensible heat flux
    add_diagnostic_variable!(
        short_name = "shf",
        long_name = "Sensible Heat Flux",
        standard_name = "sensible_heat_flux",
        units = "W m^-2",
        comments = "Exchange of energy at the land-atmosphere interface due to temperature difference.",
        compute! = (out, Y, p, t) ->
            compute_sensible_heat_flux!(out, Y, p, t, land_model),
    )

    # Vapor flux
    add_diagnostic_variable!(
        short_name = "vflux",
        long_name = "Liquid water evaporation",
        standard_name = "vapor_flux",
        units = "m s^-1",
        comments = "Flux of water from the land surface to the atmosphere. E.g., evaporation or sublimation.",
        compute! = (out, Y, p, t) ->
            compute_vapor_flux!(out, Y, p, t, land_model),
    )

    # Surface air density
    add_diagnostic_variable!(
        short_name = "rhosfc",
        long_name = "Surface Air Density",
        standard_name = "surface_air_density",
        units = "kg m^−3",
        comments = "Density of air at the land-atmosphere interface.",
        compute! = (out, Y, p, t) ->
            compute_surface_air_density!(out, Y, p, t, land_model),
    )

    # Stored in Y

    # Soil temperature (3D) at depth
    add_diagnostic_variable!(
        short_name = "tsoil",
        long_name = "Soil temperature",
        standard_name = "soil_temperature",
        units = "K",
        comments = "Soil temperature at multiple soil depth.",
        compute! = (out, Y, p, t) ->
            compute_soil_temperature!(out, Y, p, t, land_model),
    )

    # Surbsurface water storage
    add_diagnostic_variable!(
        short_name = "wsoil",
        long_name = "subsurface Water Storage",
        standard_name = "subsurface_water_storage",
        units = "m",
        comments = "Soil water content.",
        compute! = (out, Y, p, t) ->
            compute_subsurface_water_storage!(out, Y, p, t, land_model),
    )

    # Surface water content
    add_diagnostic_variable!(
        short_name = "wsfc",
        long_name = "Surface Water Content",
        standard_name = "surface_water_content",
        units = "m",
        comments = "Water at the soil surface.",
        compute! = (out, Y, p, t) ->
            compute_surface_water_content!(out, Y, p, t, land_model),
    )

    # Surface snow water content
    add_diagnostic_variable!(
        short_name = "ssfc",
        long_name = "Snow Water Equivalent",
        standard_name = "snow_water_equivalent",
        units = "m",
        comments = "Snow at the soil surface, expressed in water equivalent.",
        compute! = (out, Y, p, t) ->
            compute_snow_water_equivalent!(out, Y, p, t, land_model),
    )

    ###### SoilCanopyModel ######

    add_diagnostic_variable!(
        short_name = "slw",
        long_name = "Soil Liquid Water",
        standard_name = "soil_liquid_water",
        units = "m^3 m^-3",
        comments = "Soil moisture, the liquid water volume per soil volume. This liquid water is located in the soil pores.",
        compute! = (out, Y, p, t) ->
            compute_soil_water_liquid!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "infil",
        long_name = "Infiltration",
        standard_name = "infiltration",
        units = "m s^-1", # double check
        comments = "The flux of liquid water into the soil.",
        compute! = (out, Y, p, t) ->
            compute_infiltration!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "scd",
        long_name = "Soil CO2 Diffusivity",
        standard_name = "soil_co2_diffusivity",
        units = "", # need to add
        comments = "The diffusivity of CO2 in the porous phase of the soil. Depends on soil texture, moisture, and temperature.",
        compute! = (out, Y, p, t) ->
            compute_soilco2_diffusivity!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "scms",
        long_name = "Soil CO2 Microbial Source",
        standard_name = "soil_co2_microbial_source",
        units = "", # check
        comments = "The production of CO2 by microbes in the soil. Vary by layers of soil depth.",
        compute! = (out, Y, p, t) ->
            compute_soilco2_source_microbe!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "gs",
        long_name = "Stomatal Conductance",
        standard_name = "stomatal_conductance",
        units = "m s^-1",
        comments = "The conductance of leaves. This depends on stomatal opening. It varies with factors such as soil moisture or atmospheric water demand.",
        compute! = (out, Y, p, t) ->
            compute_stomatal_conductance!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "mt",
        long_name = "Medlyn Term",
        standard_name = "medlyn_term",
        units = "", # check
        comments = "",
        compute! = (out, Y, p, t) ->
            compute_medlyn_term!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "trans",
        long_name = "Canopy Transpiration",
        standard_name = "canopy_transpiration",
        units = "", # check
        comments = "The water evaporated from the canopy due to leaf transpiration.",
        compute! = (out, Y, p, t) ->
            compute_canopy_transpiration!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!( # not actually computed, but read from input or atmosphere model, and stored in p...
        short_name = "rain",
        long_name = "Rainfall",
        standard_name = "rainfall",
        units = "m",
        comments = "Precipitation of liquid water.",
        compute! = (out, Y, p, t) ->
            compute_rainfall!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "an",
        long_name = "Leaf Net Photosynthesis",
        standard_name = "leaf_net_photosynthesis",
        units = "", # check
        comments = "Net photosynthesis of a leaf.",
        compute! = (out, Y, p, t) ->
            compute_photosynthesis_net_leaf!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "gpp",
        long_name = "Gross Primary Productivity",
        standard_name = "gross_primary_productivity",
        units = "",
        comments = "Net photosynthesis of the canopy.",
        compute! = (out, Y, p, t) ->
            compute_photosynthesis_net_canopy!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "rd",
        long_name = "Leaf Respiration",
        standard_name = "leaf_dark_respiration",
        units = "",
        comments = "Leaf respiration, called dark respiration because usually measured in the abscence of radiation.",
        compute! = (out, Y, p, t) ->
            compute_respiration_leaf!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "vcmax25",
        long_name = "Vcmax25",
        standard_name = "vcmax25",
        units = "",
        comments = "The parameter vcmax at 25 degree celsius. Important for the Farquhar model of leaf photosynthesis.",
        compute! = (out, Y, p, t) -> compute_vcmax25!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "par",
        long_name = "Photosynthetically Active Radiation",
        standard_name = "photosynthetically_active_radiation",
        units = "",
        comments = "The subset of total radiation that activates photosynthesis reaching the canopy.",
        compute! = (out, Y, p, t) ->
            compute_photosynthetically_active_radiation!(
                out,
                Y,
                p,
                t,
                land_model,
            ),
    )

    add_diagnostic_variable!(
        short_name = "apar",
        long_name = "Absorbed Photosynthetically Active Radiation",
        standard_name = "absorbed_photosynthetically_active_radiation",
        units = "",
        comments = "The amount of photosynthetically active radiation absorbed by the leaf. The rest if reflected or transmitted.",
        compute! = (out, Y, p, t) ->
            compute_photosynthetically_active_radiation_absorbed!(
                out,
                Y,
                p,
                t,
                land_model,
            ),
    )

    add_diagnostic_variable!(
        short_name = "rpar",
        long_name = "Reflected Photosynthetically Active Radiation",
        standard_name = "reflected_photosynthetically_active_radiation",
        units = "",
        comments = "The amount of photosynthetically active radiation reflected by leaves.",
        compute! = (out, Y, p, t) ->
            compute_photosynthetically_active_radiation_reflected!(
                out,
                Y,
                p,
                t,
                land_model,
            ),
    )

    add_diagnostic_variable!(
        short_name = "tpar",
        long_name = "Transmitted Photosynthetically Active Radiation",
        standard_name = "transmitted_photosynthetically_active_radiation",
        units = "",
        comments = "The amount of photosynthetically active radiation transmitted by leaves.",
        compute! = (out, Y, p, t) ->
            compute_photosynthetically_active_radiation_transmitted!(
                out,
                Y,
                p,
                t,
                land_model,
            ),
    )

    add_diagnostic_variable!(
        short_name = "nir",
        long_name = "Near Infrared Radiation",
        standard_name = "near_infrared_radiation",
        units = "W m^-2",
        comments = "The amount of near infrared radiation reaching the canopy.",
        compute! = (out, Y, p, t) ->
            compute_near_infrared_radiation!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "anir",
        long_name = "Absorbed Near Infrared Radiation",
        standard_name = "absorbed_near_infrared_radiation",
        units = "W m^-2",
        comments = "The amount of near infrared radiation reaching the canopy.",
        compute! = (out, Y, p, t) ->
            compute_near_infrared_radiation_absorbed!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "rnir",
        long_name = "Reflected Near Infrared Radiation",
        standard_name = "reflected_near_infrared_radiation",
        units = "W m^-2",
        comments = "The amount of near infrared radiation reaching the canopy.",
        compute! = (out, Y, p, t) ->
            compute_near_infrared_radiation_reflected!(
                out,
                Y,
                p,
                t,
                land_model,
            ),
    )

    add_diagnostic_variable!(
        short_name = "tnir",
        long_name = "Transmitted Near Infrared Radiation",
        standard_name = "transmitted_near_infrared_radiation",
        units = "W m^-2",
        comments = "The amount of near infrared radiation reaching the canopy.",
        compute! = (out, Y, p, t) ->
            compute_near_infrared_radiation_transmitted!(
                out,
                Y,
                p,
                t,
                land_model,
            ),
    )

    add_diagnostic_variable!(
        short_name = "swn",
        long_name = "Net Shortwave Radiation",
        standard_name = "net_shortwave_radiation",
        units = "W m^-2",
        comments = "The net (in minus out) radiation at the surface.",
        compute! = (out, Y, p, t) ->
            compute_radiation_shortwave_net!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "lwn",
        long_name = "Net Longwave Radiation",
        standard_name = "net_longwave_radiation",
        units = "W m^-2",
        comments = "The net (in minus out) radiation at the surface.",
        compute! = (out, Y, p, t) ->
            compute_radiation_longwave_net!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "ra",
        long_name = "Autotrophic Respiration",
        standard_name = "autotrophic_respiration",
        units = "mol m^-2 s^-1",
        comments = "Canopy autotrophic respiration, the sum of leaves, stems and roots respiration.",
        compute! = (out, Y, p, t) ->
            compute_autotrophic_respiration!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "soilco2",
        long_name = "Soil CO2 concentration",
        standard_name = "soil_co2",
        units = "",
        comments = "Concentration of CO2 in the air of soil pores.",
        compute! = (out, Y, p, t) -> compute_soilco2!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "soilrn",
        long_name = "Soil Net Radiation",
        standard_name = "soil_net_radiation",
        units = "W m^-2",
        comments = "Net radiation at the soil surface.",
        compute! = (out, Y, p, t) ->
            compute_soil_net_radiation!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "soillhf",
        long_name = "Soil Latent Heat Flux",
        standard_name = "soil_Latent_Heat_Flux",
        units = "W m^-2",
        comments = "Soil evaporation.",
        compute! = (out, Y, p, t) ->
            compute_soil_latent_heat_flux!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "soilshf",
        long_name = "Soil Sensible Heat Flux",
        standard_name = "soil_sensible_Heat_Flux",
        units = "W m^-2",
        comments = "Soil sensible heat flux.",
        compute! = (out, Y, p, t) ->
            compute_soil_sensible_heat_flux!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "soilrae",
        long_name = "Soil Aerodynamic Resistance",
        standard_name = "soil_aerodynamic_resistance",
        units = "",
        comments = "Soil aerodynamic resistance.",
        compute! = (out, Y, p, t) ->
            compute_soil_aerodynamic_resistance!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "hr",
        long_name = "Heterotrophic Respiration",
        standard_name = "heterotrophic_respiration",
        units = "mol m^-2 s^-1",
        comments = "CO2 efflux at the soil surface due to microbial decomposition of soil organic matter.",
        compute! = (out, Y, p, t) ->
            compute_heterotrophic_respiration!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "shc",
        long_name = "Soil Hydraulic Conductivity",
        standard_name = "soil_hydraulic_conductivity",
        units = "",
        comments = "Soil hydraulic conductivity.",
        compute! = (out, Y, p, t) ->
            compute_soil_hydraulic_conductivity!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "stc",
        long_name = "Soil Thermal Conductivity",
        standard_name = "soil_thermal_conductivity",
        units = "",
        comments = "Soil thermal conductivity.",
        compute! = (out, Y, p, t) ->
            compute_soil_thermal_conductivity!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "swp",
        long_name = "Soil Water Potential",
        standard_name = "soil_water_potential",
        units = "",
        comments = "Soil water potential.",
        compute! = (out, Y, p, t) ->
            compute_soil_water_potential!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "sza",
        long_name = "Solar Zenith Angle",
        standard_name = "solar_zenith_angle",
        units = "",
        comments = "Solar zenith angle.",
        compute! = (out, Y, p, t) ->
            compute_solar_zenith_angle!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "msf",
        long_name = "Moisture Stress Factor",
        standard_name = "moisture_stress_factor",
        units = "",
        comments = "Sensitivity of plants conductance to soil water content.",
        compute! = (out, Y, p, t) ->
            compute_moisture_stress_factor!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "cwp",
        long_name = "Canopy Water Potential",
        standard_name = "canopy_water_potential",
        units = "",
        comments = "The water potential of the canopy.",
        compute! = (out, Y, p, t) ->
            compute_canopy_water_potential!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "fa",
        long_name = "Cross Section",
        standard_name = "cross_section",
        units = "",
        comments = "The area of stem relative to ground area.", #??
        compute! = (out, Y, p, t) ->
            compute_cross_section!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "far",
        long_name = "Root Cross Section",
        standard_name = "Cross Section",
        units = "",
        comments = "The area of roots relative to ground area.", #??
        compute! = (out, Y, p, t) ->
            compute_cross_section_roots!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "ai",
        long_name = "Area Index",
        standard_name = "area_index",
        units = "",
        comments = "The area index of leaves.",
        compute! = (out, Y, p, t) ->
            compute_area_index!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "clhf",
        long_name = "Canopy Latent Heat Flux",
        standard_name = "canopy_latent_heat_flux",
        units = "",
        comments = "Canopy evaporation.", #?? of steam, leaves, roots?
        compute! = (out, Y, p, t) ->
            compute_canopy_latent_heat_flux!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "cshf",
        long_name = "Canopy Sensible Heat Flux",
        standard_name = "canopy_sensible_heat_flux",
        units = "",
        comments = "Canopy sensible heat flux.", #?? of steam, leaves, roots?
        compute! = (out, Y, p, t) ->
            compute_canopy_sensible_heat_flux!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "crae",
        long_name = "Canopy Aerodynamic Resistance",
        standard_name = "canopy_aerodynamic_resistance",
        units = "",
        comments = "Canopy aerodynamic_resistance.", #?? of steam, leaves, roots?
        compute! = (out, Y, p, t) ->
            compute_canopy_aerodynamic_resistance!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "ct",
        long_name = "Canopy Temperature",
        standard_name = "canopy_temperature",
        units = "K",
        comments = "Canopy temperature.", #?? of steam, leaves, roots?
        compute! = (out, Y, p, t) ->
            compute_canopy_temperature!(out, Y, p, t, land_model),
    )

    add_diagnostic_variable!(
        short_name = "si",
        long_name = "Soil Ice",
        standard_name = "soil_ice",
        units = "m^3 m^-3",
        comments = "soil ice.",
        compute! = (out, Y, p, t) ->
            compute_soil_ice!(out, Y, p, t, land_model),
    )
end

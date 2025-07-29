using Dates
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput

using ClimaLand
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Canopy
import ClimaLand

"""
     global_land_model(FT,
                       context,
                       scalar_soil_params,
                       scalar_canopy_params,
                       scalar_snow_params,
                       earth_param_set;
                       context = nothing,
                       domain = ClimaLand.Domains.global_domain(FT; context = context),
                       forcing = ClimaLand.prescribed_forcing_era5(ClimaLand.Artifacts.era5_land_forcing_data2008_path(; context),
                                                                   domain.space.surface,
                                                                   DateTime(2008),
                                                                   earth_param_set,
                                                                   FT),
                       LAI = ClimaLand.prescribed_lai_modis(ClimaLand.Artifacts.modis_lai_single_year_path(DateTime(2008); context),
                                                            domain.space.surface,
                                                            DateTime(2008))
                           )

An helper function for creating a land model corresponding to a global simulation of the snow, soil, and canopy models.

While not explicitly an outer constructor, this creates and returns `LandModel` struct. This is meant as a helper for setting up a standard
global land model easily, and is useful for testing. Note that the user can construct any land model they
wish by using the default (inner) constructor method for `LandModel`, or using the alternate outer constructor method defined in src/integrated/land.jl.

Over time, all scalar parameters will be moved to ClimaParameters, so that only a single parameter set `earth_param_set` is passed.
"""
function global_land_model(
    FT,
    scalar_soil_params,
    scalar_canopy_params,
    scalar_snow_params,
    earth_param_set;
    context = nothing,
    domain = ClimaLand.Domains.global_domain(FT; context = context),
    forcing = ClimaLand.prescribed_forcing_era5(
        ClimaLand.Artifacts.era5_land_forcing_data2008_path(; context),
        domain.space.surface,
        DateTime(2008),
        earth_param_set,
        FT,
    ),
    LAI = ClimaLand.prescribed_lai_modis(
        ClimaLand.Artifacts.modis_lai_single_year_path(;
            context = nothing,
            year = 2008,
        ),
        domain.space.surface,
        DateTime(2008),
    ),
)
    # Unpack forcing
    (atmos, radiation) = forcing

    # Unpack scalar parameters
    (; α_snow, Δt) = scalar_snow_params
    (; f_over, R_sb) = scalar_soil_params
    (;
        ac_canopy,
        K_sat_plant,
        a,
        ψ63,
        Weibull_param,
        plant_ν,
        plant_S_s,
        h_leaf,
    ) = scalar_canopy_params

    # Construct spatially varying parameters.
    surface_space = domain.space.surface
    subsurface_space = domain.space.subsurface

    (; ν_ss_om, ν_ss_quartz, ν_ss_gravel) =
        ClimaLand.Soil.soil_composition_parameters(subsurface_space, FT)
    (; ν, hydrology_cm, K_sat, θ_r) =
        ClimaLand.Soil.soil_vangenuchten_parameters(subsurface_space, FT)
    soil_albedo = Soil.CLMTwoBandSoilAlbedo{FT}(;
        ClimaLand.Soil.clm_soil_albedo_parameters(surface_space)...,
    )
    S_s = ClimaCore.Fields.zeros(subsurface_space) .+ FT(1e-3)
    soil_params = Soil.EnergyHydrologyParameters(
        FT;
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm,
        K_sat,
        S_s,
        θ_r,
        albedo = soil_albedo,
    )
    f_over = FT(3.28) # 1/m
    R_sb = FT(1.484e-4 / 1000) # m/s
    runoff_model = ClimaLand.Soil.Runoff.TOPMODELRunoff{FT}(;
        f_over = f_over,
        f_max = ClimaLand.Soil.topmodel_fmax(surface_space, FT),
        R_sb = R_sb,
    )

    # Spatially varying canopy parameters from CLM
    g1 = ClimaLand.Canopy.clm_medlyn_g1(surface_space)
    rooting_depth = ClimaLand.Canopy.clm_rooting_depth(surface_space)
    (; is_c3, Vcmax25) =
        ClimaLand.Canopy.clm_photosynthesis_parameters(surface_space)
    (; Ω, G_Function, α_PAR_leaf, τ_PAR_leaf, α_NIR_leaf, τ_NIR_leaf) =
        ClimaLand.Canopy.clm_canopy_radiation_parameters(surface_space)

    conductivity_model =
        Canopy.PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
    retention_model = Canopy.PlantHydraulics.LinearRetentionCurve{FT}(a)
    SAI = FT(0.0) # m2/m2
    RAI = FT(1.0)
    n_stem = 0
    n_leaf = 1
    h_stem = FT(0.0)
    zmax = FT(0.0)

    h_canopy = h_stem + h_leaf
    compartment_midpoints =
        n_stem > 0 ? [h_stem / 2, h_stem + h_leaf / 2] : [h_leaf / 2]
    compartment_surfaces =
        n_stem > 0 ? [zmax, h_stem, h_canopy] : [zmax, h_leaf]

    z0_m = FT(0.13) * h_canopy
    z0_b = FT(0.1) * z0_m

    soil_args = (domain = domain, parameters = soil_params)
    soil_model_type = Soil.EnergyHydrology{FT}

    # Soil microbes model
    soilco2_type = Soil.Biogeochemistry.SoilCO2Model{FT}
    soilco2_ps = Soil.Biogeochemistry.SoilCO2ModelParameters(FT)
    Csom = ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))
    soilco2_args = (; domain = domain, parameters = soilco2_ps)

    # Now we set up the canopy model, which we set up by component:
    # Component Types
    canopy_component_types = (;
        autotrophic_respiration = Canopy.AutotrophicRespirationModel{FT},
        radiative_transfer = Canopy.TwoStreamModel{FT},
        photosynthesis = Canopy.FarquharModel{FT},
        conductance = Canopy.MedlynConductanceModel{FT},
        hydraulics = Canopy.PlantHydraulicsModel{FT},
        energy = Canopy.BigLeafEnergyModel{FT},
    )
    # Individual Component arguments
    # Set up autotrophic respiration
    autotrophic_respiration_args =
        (; parameters = Canopy.AutotrophicRespirationParameters(FT))
    # Set up radiative transfer
    radiative_transfer_args = (;
        parameters = Canopy.TwoStreamParameters(
            FT;
            Ω,
            α_PAR_leaf,
            τ_PAR_leaf,
            α_NIR_leaf,
            τ_NIR_leaf,
            G_Function,
        )
    )
    # Set up conductance
    conductance_args =
        (; parameters = Canopy.MedlynConductanceParameters(FT; g1))
    # Set up photosynthesis
    photosynthesis_args =
        (; parameters = Canopy.FarquharParameters(FT, is_c3; Vcmax25 = Vcmax25))
    # Set up plant hydraulics
    ai_parameterization = Canopy.PrescribedSiteAreaIndex{FT}(LAI, SAI, RAI)

    plant_hydraulics_ps = Canopy.PlantHydraulics.PlantHydraulicsParameters(;
        ai_parameterization = ai_parameterization,
        ν = plant_ν,
        S_s = plant_S_s,
        rooting_depth = rooting_depth,
        conductivity_model = conductivity_model,
        retention_model = retention_model,
    )
    plant_hydraulics_args = (
        parameters = plant_hydraulics_ps,
        n_stem = n_stem,
        n_leaf = n_leaf,
        compartment_midpoints = compartment_midpoints,
        compartment_surfaces = compartment_surfaces,
    )

    energy_args = (parameters = Canopy.BigLeafEnergyParameters{FT}(ac_canopy),)

    # Canopy component args
    canopy_component_args = (;
        autotrophic_respiration = autotrophic_respiration_args,
        radiative_transfer = radiative_transfer_args,
        photosynthesis = photosynthesis_args,
        conductance = conductance_args,
        hydraulics = plant_hydraulics_args,
        energy = energy_args,
    )

    # Other info needed
    shared_params = Canopy.SharedCanopyParameters{FT, typeof(earth_param_set)}(
        z0_m,
        z0_b,
        earth_param_set,
    )

    canopy_model_args = (;
        parameters = shared_params,
        domain = ClimaLand.obtain_surface_domain(domain),
    )

    # Snow model
    snow_parameters = SnowParameters{FT}(
        Δt;
        earth_param_set = earth_param_set,
        α_snow = α_snow,
    )
    snow_args = (;
        parameters = snow_parameters,
        domain = ClimaLand.obtain_surface_domain(domain),
    )
    snow_model_type = Snow.SnowModel

    land_input = (
        atmos = atmos,
        radiation = radiation,
        runoff = runoff_model,
        soil_organic_carbon = Csom,
    )
    return LandModel{FT}(;
        soilco2_type = soilco2_type,
        soilco2_args = soilco2_args,
        land_args = land_input,
        soil_model_type = soil_model_type,
        soil_args = soil_args,
        canopy_component_types = canopy_component_types,
        canopy_component_args = canopy_component_args,
        canopy_model_args = canopy_model_args,
        snow_args = snow_args,
        snow_model_type = snow_model_type,
    )
end

"""
    check_ocean_values_p(p, binary_mask; val = 0.0)

This function tests that every field stored in `p` has all of
its values (where binary_mask == 1) equal to  `val`. Note that
this is meant to be used with the full land model (canopy,
snow, soil, soilco2).

Useful for checking if land model functions are updating the
values over the ocean.
"""
function check_ocean_values_p(p, binary_mask; val = 0.0)
    properties = [
        p.drivers,
        p.soil,
        p.soilco2,
        p.snow,
        p.canopy.energy,
        p.canopy.hydraulics,
        p.canopy.radiative_transfer,
        p.canopy.photosynthesis,
        p.canopy.sif,
        p.canopy.turbulent_fluxes,
        p.canopy.autotrophic_respiration,
        p.canopy.conductance,
    ]
    for property in properties
        for var in propertynames(property)
            field_values = Array(parent(getproperty(property, var)))
            if length(size(field_values)) == 5 # 3d var
                @test extrema(field_values[:, 1, 1, :, Array(binary_mask)]) ==
                      (val, val)
            else
                @test extrema(field_values[1, 1, :, Array(binary_mask)]) ==
                      (val, val)
            end
        end
    end

    field_pn_p = [
        pn for pn in propertynames(p) if pn != :soil &&
        pn != :canopy &&
        pn != :snow &&
        pn != :soilco2 &&
        pn != :drivers &&
        ~occursin("dss", String(pn))
    ]

    for var in field_pn_p
        field_values = Array(parent(getproperty(p, var)))
        if length(size(field_values)) == 5 # 3d var
            @test extrema(field_values[:, 1, 1, :, Array(binary_mask)]) ==
                  (val, val)
        else
            @test extrema(field_values[1, 1, :, Array(binary_mask)]) ==
                  (val, val)

        end
    end
end

"""
    check_ocean_values_Y(Y, binary_mask; val = 0.0)

This function tests that every field stored in `Y` has all of
its values (where binary_mask == 1) equal to  `val`. Note that
this is meant to be used with the full land model (canopy,
snow, soil, soilco2).

Useful for checking if land model functions are updating the
values over the ocean.
"""
function check_ocean_values_Y(Y, binary_mask; val = 0.0)
    @test extrema(Array(parent(Y.soil.ϑ_l))[:, 1, 1, 1, Array(binary_mask)]) ==
          (val, val)
    @test extrema(Array(parent(Y.soil.θ_i))[:, 1, 1, 1, Array(binary_mask)]) ==
          (val, val)
    @test extrema(
        Array(parent(Y.soil.ρe_int))[:, 1, 1, 1, Array(binary_mask)],
    ) == (val, val)
    @test extrema(Array(parent(Y.snow.U))[1, 1, 1, Array(binary_mask)]) ==
          (val, val)
    @test extrema(Array(parent(Y.snow.S))[1, 1, 1, Array(binary_mask)]) ==
          (val, val)
    @test extrema(Array(parent(Y.snow.S_l))[1, 1, 1, Array(binary_mask)]) ==
          (val, val)
    @test extrema(
        Array(parent(Y.canopy.energy.T))[1, 1, 1, Array(binary_mask)],
    ) == (val, val)
    @test extrema(
        Array(parent(Y.canopy.hydraulics.ϑ_l.:1))[1, 1, 1, Array(binary_mask)],
    ) == (val, val)
end

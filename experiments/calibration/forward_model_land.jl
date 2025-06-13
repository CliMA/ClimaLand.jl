ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON" # force climacomms to use singleton context and not MPI

import SciMLBase
import TOML
import ClimaComms
ClimaComms.@import_required_backends
import ClimaTimeSteppers as CTS
using ClimaCore
import ClimaUtilities.OnlineLogging: WallTimeInfo, report_walltime
import ClimaUtilities.TimeManager: ITime, date

import ClimaDiagnostics
import ClimaUtilities

import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP

using ClimaLand
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Canopy
import ClimaLand.Parameters as LP

import ClimaCalibrate: forward_model, parameter_path, path_to_ensemble_member
import ClimaCalibrate as CAL

using Statistics
using Dates
import NCDatasets

FT = Float64
const time_interpolation_method = LinearInterpolation()
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"

function setup_model(
    FT,
    start_date,
    stop_date,
    Δt,
    domain,
    earth_param_set,
    calibration_params,
)
    (; α_0, k, Δα) = calibration_params

    surface_space = domain.space.surface
    subsurface_space = domain.space.subsurface
    # Forcing data
    era5_ncdata_path =
        ClimaLand.Artifacts.find_era5_year_paths(start_date, stop_date; context)
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        era5_ncdata_path,
        surface_space,
        start_date,
        earth_param_set,
        FT;
        max_wind_speed = 25.0,
        time_interpolation_method,
    )

    (; ν_ss_om, ν_ss_quartz, ν_ss_gravel) =
        ClimaLand.Soil.soil_composition_parameters(subsurface_space, FT)
    (; ν, hydrology_cm, K_sat, θ_r) =
        ClimaLand.Soil.soil_vangenuchten_parameters(subsurface_space, FT)
    soil_albedo = Soil.CLMTwoBandSoilAlbedo{FT}(;
        ClimaLand.Soil.clm_soil_albedo_parameters(surface_space)...,)
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

    # Energy Balance model
    ac_canopy = FT(2.5e3)
    # Plant Hydraulics and general plant parameters
    SAI = FT(0.0) # m2/m2
    f_root_to_shoot = FT(3.5)
    RAI = FT(1.0)
    K_sat_plant = FT(7e-8) # m/s 
    ψ63 = FT(-4 / 0.0098) # / MPa to m
    Weibull_param = FT(4) # unitless
    a = FT(0.2 * 0.0098) # 1/m
    conductivity_model =
        Canopy.PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
    retention_model = Canopy.PlantHydraulics.LinearRetentionCurve{FT}(a)
    plant_ν = FT(1.44e-4)
    plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
    n_stem = 0
    n_leaf = 1
    h_stem = FT(0.0)
    h_leaf = FT(1.0)
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
    modis_lai_ncdata_path = ClimaLand.Artifacts.find_modis_year_paths(
        start_date,
        stop_date;
        context,
    )

    LAIfunction = ClimaLand.prescribed_lai_modis(
        modis_lai_ncdata_path,
        surface_space,
        start_date;
        time_interpolation_method,
    )
    ai_parameterization =
        Canopy.PrescribedSiteAreaIndex{FT}(LAIfunction, SAI, RAI)

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
    α_snow = Snow.ZenithAngleAlbedoModel(FT(α_0), FT(Δα), FT(k))

    horz_degree_res =
        sum(ClimaLand.Domains.average_horizontal_resolution_degrees(domain)) / 2 # mean of resolution in latitude and longitude, in degrees
    scf = Snow.WuWuSnowCoverFractionModel(
        FT(0.08),
        FT(1.77),
        FT(1.0),
        horz_degree_res
    )
    snow_parameters = SnowParameters{FT}(
        Δt;
        earth_param_set = earth_param_set,
        α_snow = α_snow,
        scf = scf,
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
    land = LandModel{FT}(;
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
    return land
end

function CAL.forward_model(iteration, member)
    ensemble_member_path = path_to_ensemble_member(caldir, iteration, member)
    params_path = parameter_path(caldir, iteration, member)
    params = TOML.parsefile(params_path)

    sim_start = start_date + Year(iteration)
    sim_end = start_date + Year(iteration) + Year(2) + Day(1) # need to go slightly past two years to get the output written from the last month
    Δt = 450.0

    domain =
        ClimaLand.ModelSetup.global_domain(FT; comms_ctx = context, nelements)
    earth_param_set = LP.LandParameters(FT)
    p_names = collect(keys(params))
    p_values = [params[name]["value"] for name in p_names]
    calibration_params = (; zip(Symbol.(p_names), p_values)...)
    model = setup_model(
        FT,
        sim_start,
        sim_end,
        Δt,
        domain,
        earth_param_set,
        calibration_params,
    )

    diagnostics_dir = joinpath(ensemble_member_path, "global_diagnostics")
    diagdir =
        ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_dir)
    nc_writer = ClimaDiagnostics.Writers.NetCDFWriter(
        domain.space.subsurface,
        diagdir;
        start_date = sim_start,
    )
    diagnostics = ClimaLand.default_diagnostics(
        model,
        sim_start;
        output_writer = nc_writer,
    )

    simulation = ClimaLand.Simulations.LandSimulation(
        FT,
        sim_start,
        sim_end,
        Δt,
        model;
	user_callbacks = (),
        outdir = diagdir,
        diagnostics,
    )
    ClimaLand.Simulations.solve!(simulation)
    close(nc_writer)
    return nothing
end

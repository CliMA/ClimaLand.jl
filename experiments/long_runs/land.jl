# # Global run of land model

# The code sets up and runs the soil/canopy model for 6 hours on a spherical domain,
# using ERA5 data. In this simulation, we have
# turned lateral flow off because horizontal boundary conditions and the
# land/sea mask are not yet supported by ClimaCore.

# Simulation Setup
# Number of spatial elements: 101 in horizontal, 15 in vertical
# Soil depth: 50 m
# Simulation duration: 365 d
# Timestep: 900 s
# Timestepper: ARS343
# Fixed number of iterations: 1
# Jacobian update: every new timestep
# Atmos forcing update: every 3 hours
import SciMLBase
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import ClimaTimeSteppers as CTS
using ClimaCore
using ClimaUtilities.ClimaArtifacts
import Interpolations
using Insolation

import ClimaDiagnostics
import ClimaAnalysis
import ClimaAnalysis.Visualize as viz
import ClimaUtilities

import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP

using ClimaLand
using ClimaLand.Soil
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP

using Statistics
using CairoMakie
import GeoMakie
using Dates
import NCDatasets

const FT = Float64;

regridder_type = :InterpolationsRegridder
context = ClimaComms.context()
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "land_longrun_$(device_suffix)"
diagnostics_outdir = joinpath(root_path, "global_diagnostics")
outdir =
    ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir)

function setup_prob(t0, tf, Δt; outdir = outdir, nelements = (101, 15))

    earth_param_set = LP.LandParameters(FT)
    radius = FT(6378.1e3)
    depth = FT(50)
    domain = ClimaLand.Domains.SphericalShell(;
        radius = radius,
        depth = depth,
        nelements = nelements,
        npolynomial = 1,
        dz_tuple = FT.((10.0, 0.5)),# top layer should ideally be only a few cm!
    )
    surface_space = domain.space.surface
    subsurface_space = domain.space.subsurface

    ref_time = DateTime(2021)
    t_start = t0
    # Forcing data
    era5_artifact_path =
        ClimaLand.Artifacts.era5_land_forcing_data2021_folder_path(; context)    # Precipitation:
    precip = TimeVaryingInput(
        joinpath(era5_artifact_path, "era5_2021_0.9x1.25_clima.nc"),
        "rf",
        surface_space;
        reference_date = ref_time,
        t_start,
        regridder_type,
        file_reader_kwargs = (; preprocess_func = (data) -> -data / 3600,),
    )

    snow_precip = TimeVaryingInput(
        joinpath(era5_artifact_path, "era5_2021_0.9x1.25.nc"),
        "sf",
        surface_space;
        reference_date = ref_time,
        t_start,
        regridder_type,
        file_reader_kwargs = (; preprocess_func = (data) -> -data / 3600,),
    )

    u_atmos = TimeVaryingInput(
        joinpath(era5_artifact_path, "era5_2021_0.9x1.25_clima.nc"),
        "ws",
        surface_space;
        reference_date = ref_time,
        t_start,
        regridder_type,
    )
    q_atmos = TimeVaryingInput(
        joinpath(era5_artifact_path, "era5_2021_0.9x1.25_clima.nc"),
        "q",
        surface_space;
        reference_date = ref_time,
        t_start,
        regridder_type,
    )
    P_atmos = TimeVaryingInput(
        joinpath(era5_artifact_path, "era5_2021_0.9x1.25.nc"),
        "sp",
        surface_space;
        reference_date = ref_time,
        t_start,
        regridder_type,
    )

    T_atmos = TimeVaryingInput(
        joinpath(era5_artifact_path, "era5_2021_0.9x1.25.nc"),
        "t2m",
        surface_space;
        reference_date = ref_time,
        t_start,
        regridder_type,
    )
    h_atmos = FT(10)

    atmos = PrescribedAtmosphere(
        precip,
        snow_precip,
        T_atmos,
        u_atmos,
        q_atmos,
        P_atmos,
        ref_time,
        h_atmos,
        earth_param_set,
    )

    # Prescribed radiation
    SW_d = TimeVaryingInput(
        joinpath(era5_artifact_path, "era5_2021_0.9x1.25.nc"),
        "ssrd",
        surface_space;
        reference_date = ref_time,
        t_start,
        regridder_type,
        file_reader_kwargs = (; preprocess_func = (data) -> data / 3600,),
    )
    LW_d = TimeVaryingInput(
        joinpath(era5_artifact_path, "era5_2021_0.9x1.25.nc"),
        "strd",
        surface_space;
        reference_date = ref_time,
        t_start,
        regridder_type,
        file_reader_kwargs = (; preprocess_func = (data) -> data / 3600,),
    )

    function zenith_angle(
        t,
        ref_time;
        latitude = ClimaCore.Fields.coordinate_field(surface_space).lat,
        longitude = ClimaCore.Fields.coordinate_field(surface_space).long,
        insol_params::Insolation.Parameters.InsolationParameters{FT} = earth_param_set.insol_params,
    ) where {FT}
        # This should be time in UTC
        current_datetime = ref_time + Dates.Second(round(t))

        # Orbital Data uses Float64, so we need to convert to our sim FT
        d, δ, η_UTC =
            FT.(
                Insolation.helper_instantaneous_zenith_angle(
                    current_datetime,
                    ref_time,
                    insol_params,
                )
            )

        Insolation.instantaneous_zenith_angle.(
            d,
            δ,
            η_UTC,
            longitude,
            latitude,
        ).:1
    end
    radiation =
        PrescribedRadiativeFluxes(FT, SW_d, LW_d, ref_time; θs = zenith_angle)

    soil_params_artifact_path =
        ClimaLand.Artifacts.soil_params_artifact_folder_path(; context)
    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    )
    soil_params_mask = SpaceVaryingInput(
        joinpath(
            soil_params_artifact_path,
            "vGalpha_map_gupta_etal2020_1.0x1.0x4.nc",
        ),
        "α",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
        file_reader_kwargs = (; preprocess_func = (data) -> data > 0,),
    )
    oceans_to_value(field, mask, value) =
        mask == 1.0 ? field : eltype(field)(value)

    vg_α = SpaceVaryingInput(
        joinpath(
            soil_params_artifact_path,
            "vGalpha_map_gupta_etal2020_1.0x1.0x4.nc",
        ),
        "α",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )
    vg_n = SpaceVaryingInput(
        joinpath(
            soil_params_artifact_path,
            "vGn_map_gupta_etal2020_1.0x1.0x4.nc",
        ),
        "n",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )
    x = parent(vg_α)
    μ = mean(log10.(x[x .> 0]))
    vg_α .= oceans_to_value.(vg_α, soil_params_mask, 10.0^μ)

    x = parent(vg_n)
    μ = mean(x[x .> 0])
    vg_n .= oceans_to_value.(vg_n, soil_params_mask, μ)

    vg_fields_to_hcm_field(α::FT, n::FT) where {FT} =
        ClimaLand.Soil.vanGenuchten{FT}(; @NamedTuple{α::FT, n::FT}((α, n))...)
    hydrology_cm = vg_fields_to_hcm_field.(vg_α, vg_n)

    θ_r = SpaceVaryingInput(
        joinpath(
            soil_params_artifact_path,
            "residual_map_gupta_etal2020_1.0x1.0x4.nc",
        ),
        "θ_r",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )

    ν = SpaceVaryingInput(
        joinpath(
            soil_params_artifact_path,
            "porosity_map_gupta_etal2020_1.0x1.0x4.nc",
        ),
        "ν",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )
    K_sat = SpaceVaryingInput(
        joinpath(
            soil_params_artifact_path,
            "ksat_map_gupta_etal2020_1.0x1.0x4.nc",
        ),
        "Ksat",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )

    x = parent(K_sat)
    μ = mean(log10.(x[x .> 0]))
    K_sat .= oceans_to_value.(K_sat, soil_params_mask, 10.0^μ)

    ν .= oceans_to_value.(ν, soil_params_mask, 1)

    θ_r .= oceans_to_value.(θ_r, soil_params_mask, 0)


    S_s =
        oceans_to_value.(
            ClimaCore.Fields.zeros(subsurface_space) .+ FT(1e-3),
            soil_params_mask,
            1,
        )
    ν_ss_om =
        oceans_to_value.(
            ClimaCore.Fields.zeros(subsurface_space) .+ FT(0.1),
            soil_params_mask,
            0,
        )
    ν_ss_quartz =
        oceans_to_value.(
            ClimaCore.Fields.zeros(subsurface_space) .+ FT(0.1),
            soil_params_mask,
            0,
        )
    ν_ss_gravel =
        oceans_to_value.(
            ClimaCore.Fields.zeros(subsurface_space) .+ FT(0.1),
            soil_params_mask,
            0,
        )
    PAR_albedo = ClimaCore.Fields.zeros(surface_space) .+ FT(0.2)
    NIR_albedo = ClimaCore.Fields.zeros(surface_space) .+ FT(0.2)
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
        PAR_albedo = PAR_albedo,
        NIR_albedo = NIR_albedo,
    )

    soil_params_mask_sfc =
        ClimaLand.Domains.top_center_to_surface(soil_params_mask)

    # Read in f_max data and land sea mask
    infile_path = ClimaLand.Artifacts.topmodel_data_path()
    f_max =
        SpaceVaryingInput(infile_path, "fmax", surface_space; regridder_type)
    mask = SpaceVaryingInput(
        infile_path,
        "landsea_mask",
        surface_space;
        regridder_type,
    )
    # Unsure how to handle two masks
    f_max = oceans_to_value.(f_max, mask, FT(0.0))
    f_max = oceans_to_value.(f_max, soil_params_mask_sfc, FT(0.0))
    f_over = FT(3.28) # 1/m
    R_sb = FT(1.484e-4 / 1000) # m/s
    runoff_model = ClimaLand.Soil.Runoff.TOPMODELRunoff{FT}(;
        f_over = f_over,
        f_max = f_max,
        R_sb = R_sb,
    )


    # TwoStreamModel parameters
    Ω = FT(0.69)
    ld = FT(0.5)
    α_PAR_leaf = FT(0.1)
    τ_PAR_leaf = FT(0.05)
    α_NIR_leaf = FT(0.45)
    τ_NIR_leaf = FT(0.25)

    # Energy Balance model
    ac_canopy = FT(2.5e4) # this will likely be 10x smaller!

    # Conductance Model
    g1 = FT(141) # Wang et al: 141 sqrt(Pa) for Medlyn model; Natan used 300.

    #Photosynthesis model
    clm_artifact_path = ClimaLand.Artifacts.clm_data_folder_path(; context)
    # vcmax is read in units of umol CO2/m^2/s and then converted to mol CO2/m^2/s
    Vcmax25 = SpaceVaryingInput(
        joinpath(clm_artifact_path, "vegetation_properties_map.nc"),
        "vcmx25",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
        file_reader_kwargs = (; preprocess_func = (data) -> data / 1_000_000,),
    )

    # Plant Hydraulics and general plant parameters
    SAI = FT(0.0) # m2/m2
    f_root_to_shoot = FT(3.5)
    RAI = FT(1.0)
    K_sat_plant = FT(5e-9) # m/s # seems much too small?
    ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value is -4 MPa
    Weibull_param = FT(4) # unitless, Holtzman's original c param value
    a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
    conductivity_model =
        Canopy.PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
    retention_model = Canopy.PlantHydraulics.LinearRetentionCurve{FT}(a)
    plant_ν = FT(1.44e-4)
    plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
    rooting_depth = FT(0.5) # from Natan
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


    soilco2_ps = Soil.Biogeochemistry.SoilCO2ModelParameters(FT)

    soil_args = (domain = domain, parameters = soil_params)
    soil_model_type = Soil.EnergyHydrology{FT}

    # Soil microbes model
    soilco2_type = Soil.Biogeochemistry.SoilCO2Model{FT}

    # soil microbes args
    Csom = ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))

    # Set the soil CO2 BC to being atmospheric CO2
    soilco2_top_bc = Soil.Biogeochemistry.AtmosCO2StateBC()
    soilco2_bot_bc = Soil.Biogeochemistry.SoilCO2FluxBC((p, t) -> 0.0) # no flux
    soilco2_sources = (Soil.Biogeochemistry.MicrobeProduction{FT}(),)

    soilco2_boundary_conditions =
        (; top = soilco2_top_bc, bottom = soilco2_bot_bc)

    soilco2_args = (;
        boundary_conditions = soilco2_boundary_conditions,
        sources = soilco2_sources,
        domain = domain,
        parameters = soilco2_ps,
    )

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
        )
    )
    # Set up conductance
    conductance_args =
        (; parameters = Canopy.MedlynConductanceParameters(FT; g1))
    # Set up photosynthesis
    photosynthesis_args = (;
        parameters = Canopy.FarquharParameters(
            FT,
            Canopy.C3();
            Vcmax25 = Vcmax25,
        )
    )
    # Set up plant hydraulics

    # Note that we clip all values of LAI below 0.05 to zero.
    # This is because we currently run into issues when LAI is
    # of order eps(FT) in the SW radiation code.
    # Please see Issue #644
    # or PR #645 for details.
    # For now, this clipping is similar to what CLM does.
    LAIfunction = TimeVaryingInput(
        joinpath(era5_artifact_path, "era5_lai_2021_0.9x1.25_clima.nc"),
        "lai",
        surface_space;
        reference_date = ref_time,
        t_start,
        regridder_type,
        file_reader_kwargs = (;
            preprocess_func = (data) -> data > 0.05 ? data : 0.0,
        ),
    )
    ai_parameterization =
        Canopy.PrescribedSiteAreaIndex{FT}(LAIfunction, SAI, RAI)

    function root_distribution(z::T; rooting_depth = rooting_depth) where {T}
        return T(1.0 / rooting_depth) * exp(z / T(rooting_depth)) # 1/m
    end

    plant_hydraulics_ps = Canopy.PlantHydraulics.PlantHydraulicsParameters(;
        ai_parameterization = ai_parameterization,
        ν = plant_ν,
        S_s = plant_S_s,
        root_distribution = root_distribution,
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

    # Integrated plant hydraulics and soil model
    land_input = (
        atmos = atmos,
        radiation = radiation,
        runoff = runoff_model,
        soil_organic_carbon = Csom,
    )
    land = SoilCanopyModel{FT}(;
        soilco2_type = soilco2_type,
        soilco2_args = soilco2_args,
        land_args = land_input,
        soil_model_type = soil_model_type,
        soil_args = soil_args,
        canopy_component_types = canopy_component_types,
        canopy_component_args = canopy_component_args,
        canopy_model_args = canopy_model_args,
    )

    Y, p, cds = initialize(land)

    init_soil(ν, θ_r) = θ_r + (ν - θ_r) / 2
    Y.soil.ϑ_l .= init_soil.(ν, θ_r)
    Y.soil.θ_i .= FT(0.0)
    T = FT(276.85)
    ρc_s =
        Soil.volumetric_heat_capacity.(
            Y.soil.ϑ_l,
            Y.soil.θ_i,
            soil_params.ρc_ds,
            soil_params.earth_param_set,
        )
    Y.soil.ρe_int .=
        Soil.volumetric_internal_energy.(
            Y.soil.θ_i,
            ρc_s,
            T,
            soil_params.earth_param_set,
        )
    Y.soilco2.C .= FT(0.000412) # set to atmospheric co2, mol co2 per mol air
    Y.canopy.hydraulics.ϑ_l.:1 .= plant_ν
    evaluate!(Y.canopy.energy.T, atmos.T, t0)

    set_initial_cache! = make_set_initial_cache(land)
    exp_tendency! = make_exp_tendency(land)
    imp_tendency! = ClimaLand.make_imp_tendency(land)
    jacobian! = ClimaLand.make_jacobian(land)
    set_initial_cache!(p, Y, t0)

    # set up jacobian info
    jac_kwargs =
        (; jac_prototype = ImplicitEquationJacobian(Y), Wfact = jacobian!)

    prob = SciMLBase.ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = exp_tendency!,
            T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
            dss! = ClimaLand.dss!,
        ),
        Y,
        (t0, tf),
        p,
    )

    updateat = Array(t0:(3600 * 3):tf)
    drivers = ClimaLand.get_drivers(land)
    updatefunc = ClimaLand.make_update_drivers(drivers)

    # ClimaDiagnostics

    nc_writer = ClimaDiagnostics.Writers.NetCDFWriter(subsurface_space, outdir)

    diags = ClimaLand.default_diagnostics(
        land,
        t0,
        ref_time;
        output_writer = nc_writer,
        output_vars = :long,
    )

    diagnostic_handler =
        ClimaDiagnostics.DiagnosticsHandler(diags, Y, p, t0; dt = Δt)

    diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler)

    driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
    return prob, SciMLBase.CallbackSet(driver_cb, diag_cb)
end

function setup_and_solve_problem(; greet = false)

    t0 = 0.0
    tf = 60 * 60.0 * 24 * 7 # keep short until it runs! * 365
    Δt = 900.0
    nelements = (101, 15)
    if greet
        @info "Run: Global Soil-Canopy Model"
        @info "Resolution: $nelements"
        @info "Timestep: $Δt s"
        @info "Duration: $(tf - t0) s"
    end

    prob, cb = setup_prob(t0, tf, Δt; nelements)

    # Define timestepper and ODE algorithm
    stepper = CTS.ARS343()
    ode_algo = CTS.IMEXAlgorithm(
        stepper,
        CTS.NewtonsMethod(
            max_iters = 1,
            update_j = CTS.UpdateEvery(CTS.NewTimeStep),
        ),
    )
    SciMLBase.solve(prob, ode_algo; dt = Δt, callback = cb, adaptive = false)
    return nothing
end

setup_and_solve_problem(; greet = true);
# read in diagnostics and make some plots!
#### ClimaAnalysis ####
simdir = ClimaAnalysis.SimDir(outdir)
short_names = ["gpp", "swc", "si", "sie"]
for short_name in short_names
    var = get(simdir; short_name)
    times = ClimaAnalysis.times(var)
    for t in times
        fig = CairoMakie.Figure(size = (800, 600))
        kwargs = ClimaAnalysis.has_altitude(var) ? Dict(:z => 1) : Dict()
        viz.heatmap2D_on_globe!(
            fig,
            ClimaAnalysis.slice(var, time = t; kwargs...),
        )
        CairoMakie.save(joinpath(root_path, "$(short_name)_$t.png"), fig)
    end
end

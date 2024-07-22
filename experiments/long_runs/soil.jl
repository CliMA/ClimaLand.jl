# # Global run of soil model

# The code sets up and runs the soil model for on a spherical domain,
# using ERA5 data. In this simulation, we have
# turned lateral flow off because horizontal boundary conditions and the
# land/sea mask are not yet supported by ClimaCore.

# Simulation Setup
# Number of spatial elements: 101 1in horizontal, 15 in vertical
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

using ClimaDiagnostics
using ClimaAnalysis
import ClimaAnalysis.Visualize as viz
using ClimaUtilities

import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP

using ClimaLand
using ClimaLand.Soil
import ClimaLand
import ClimaLand.Parameters as LP

using Statistics
using CairoMakie
using Dates
import NCDatasets

const FT = Float64;

regridder_type = :InterpolationsRegridder
context = ClimaComms.context()
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "soil_longrun_$(device_suffix)"
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
        dz_tuple = FT.((10.0, 0.05)),
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

    soil_args = (domain = domain, parameters = soil_params)
    soil_model_type = Soil.EnergyHydrology{FT}
    sources = (Soil.PhaseChange{FT}(FT(3600)),)# sublimation and subsurface runoff are added automatically
    top_bc = ClimaLand.Soil.AtmosDrivenFluxBC(atmos, radiation, runoff_model)
    zero_flux = Soil.HeatFluxBC((p, t) -> 0.0)
    boundary_conditions = (;
        top = top_bc,
        bottom = Soil.WaterHeatBC(;
            water = Soil.FreeDrainage(),
            heat = zero_flux,
        ),
    )
    soil = soil_model_type(;
        boundary_conditions = boundary_conditions,
        sources = sources,
        soil_args...,
    )


    Y, p, cds = initialize(soil)

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

    set_initial_cache! = make_set_initial_cache(soil)
    exp_tendency! = make_exp_tendency(soil)
    imp_tendency! = ClimaLand.make_imp_tendency(soil)
    jacobian! = ClimaLand.make_jacobian(soil)
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
    drivers = ClimaLand.get_drivers(soil)
    updatefunc = ClimaLand.make_update_drivers(drivers)

    # ClimaDiagnostics

    nc_writer = ClimaDiagnostics.Writers.NetCDFWriter(subsurface_space, outdir)

    diags =
        ClimaLand.CLD.default_diagnostics(soil, t0; output_writer = nc_writer)

    diagnostic_handler =
        ClimaDiagnostics.DiagnosticsHandler(diags, Y, p, t0; dt = Δt)

    diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler)

    driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
    return prob, SciMLBase.CallbackSet(driver_cb, diag_cb)
end

function setup_and_solve_problem(; greet = false)

    t0 = 0.0
    tf = 60 * 60.0 * 24 * 60 # keep short until it runs! * 365
    Δt = 900.0
    nelements = (101, 15)
    if greet
        @info "Run: Global Soil Model"
        @info "Resolution: $nelements"
        @info "Timestep: $Δt s"
        @info "Duration: $(tf - t0) s"
    end

    prob, cb = setup_prob(t0, tf, Δt; nelements)

    # Define timestepper and ODE algorithm
    stepper = CTS.ARS111()
    ode_algo = CTS.IMEXAlgorithm(
        stepper,
        CTS.NewtonsMethod(
            max_iters = 5,
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
short_names_2D = ["slw", "si", "tsoil"]
times = 0.0:(60.0 * 60.0 * 24 * 20):(60.0 * 60.0 * 24 * 60)
var_limits = [(0, 1), (0, 1), (270, 320)]
for t in times
    for (short_name, limits) in zip(short_names_2D, var_limits)
        var = get(simdir; short_name)
        fig = CairoMakie.Figure(size = (800, 600))
        more_kwargs = Dict(:plot => Dict(:colorrange => limits))
        kwargs = Dict(:time => t)
        viz.plot!(fig, var; more_kwargs = more_kwargs, time = t)
        CairoMakie.save(joinpath(root_path, "$short_name $t.png"), fig)
    end
end

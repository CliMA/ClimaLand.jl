# # Global run of land model

# The code sets up and runs ClimaLand v1, which
# includes soil, canopy, and snow, on a spherical domain,
# using ERA5 data as forcing. In this simulation, we have
# turned lateral flow off because horizontal boundary conditions and the
# land/sea mask are not yet supported by ClimaCore.

# Simulation Setup
# Number of spatial elements: 101 in horizontal, 15 in vertical
# Soil depth: 50 m
# Simulation duration: 730 d
# Timestep: 450 s
# Timestepper: ARS111
# Fixed number of iterations: 3
# Jacobian update: every new Newton iteration
# Atmos forcing update: every 3 hours
import SciMLBase
import ClimaComms
ClimaComms.@import_required_backends
import ClimaTimeSteppers as CTS
import ClimaCore
@show pkgversion(ClimaCore)
using ClimaUtilities.ClimaArtifacts
import ClimaUtilities.OnlineLogging: WallTimeInfo, report_walltime
import ClimaUtilities.TimeManager: ITime, date

import ClimaDiagnostics
import ClimaAnalysis
import ClimaAnalysis.Visualize as viz
import ClimaUtilities

import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP

using ClimaLand
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP

using Statistics
using CairoMakie
import GeoMakie
using Dates
import NCDatasets

using Poppler_jll: pdfunite

const FT = Float64;
# If you want to do a very long run locally, you can enter `export
# LONGER_RUN=""` in the terminal and run this script. If you want to do a very
# long run on Buildkite manually, then make a new build and pass `LONGER_RUN=""`
# as an environment variable. In both cases, the value of `LONGER_RUN` does not
# matter.
const LONGER_RUN = haskey(ENV, "LONGER_RUN") ? true : false
# time_interpolation_method =
#     LONGER_RUN ? LinearInterpolation() : LinearInterpolation(PeriodicCalendar())
time_interpolation_method =
    LONGER_RUN ? LinearInterpolation(PeriodicCalendar()) : LinearInterpolation(PeriodicCalendar())
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "snowy_land_longrun_$(device_suffix)"
diagnostics_outdir = joinpath(root_path, "global_diagnostics")
outdir =
    ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir)

function setup_prob(
    t0,
    tf,
    Δt,
    start_date;
    outdir = outdir,
    nelements = (101, 15),
)
    earth_param_set = LP.LandParameters(FT)
    domain = ClimaLand.global_domain(FT; nelements = nelements)
    surface_space = domain.space.surface
    subsurface_space = domain.space.subsurface

    # Forcing data
    if LONGER_RUN
        era5_ncdata_path = ClimaLand.Artifacts.find_era5_year_paths(
            date(t0),
            date(tf);
            context,
        )
    else
        era5_artifact_path =
            ClimaLand.Artifacts.era5_land_forcing_data2008_folder_path(;
                context,
            )
        era5_ncdata_path = joinpath(era5_artifact_path, "era5_2008_1.0x1.0.nc")
    end
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        era5_ncdata_path,
        surface_space,
        start_date,
        earth_param_set,
        FT;
        time_interpolation_method = time_interpolation_method,
    )

    spatially_varying_soil_params =
        ClimaLand.default_spatially_varying_soil_parameters(
            subsurface_space,
            surface_space,
            FT,
        )
    (;
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm,
        K_sat,
        S_s,
        θ_r,
        PAR_albedo_dry,
        NIR_albedo_dry,
        PAR_albedo_wet,
        NIR_albedo_wet,
        f_max,
    ) = spatially_varying_soil_params
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
        PAR_albedo_dry = PAR_albedo_dry,
        NIR_albedo_dry = NIR_albedo_dry,
        PAR_albedo_wet = PAR_albedo_wet,
        NIR_albedo_wet = NIR_albedo_wet,
    )
    f_over = FT(3.28) # 1/m
    R_sb = FT(1.484e-4 / 1000) # m/s
    runoff_model = ClimaLand.Soil.Runoff.TOPMODELRunoff{FT}(;
        f_over = f_over,
        f_max = f_max,
        R_sb = R_sb,
    )

    # Spatially varying canopy parameters from CLM
    clm_parameters = ClimaLand.clm_canopy_parameters(surface_space)
    (;
        Ω,
        rooting_depth,
        is_c3,
        Vcmax25,
        g1,
        G_Function,
        α_PAR_leaf,
        τ_PAR_leaf,
        α_NIR_leaf,
        τ_NIR_leaf,
    ) = clm_parameters

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
    if LONGER_RUN
        modis_lai_ncdata_path = ClimaLand.Artifacts.find_modis_year_paths(
            date(t0),
            date(tf);
            context,
        )
    else
        modis_lai_artifact_path =
            ClimaLand.Artifacts.modis_lai_forcing_data_path(; context)
        modis_lai_ncdata_path =
            joinpath(modis_lai_artifact_path, "Yuan_et_al_2008_1x1.nc")
    end
    LAIfunction = ClimaLand.prescribed_lai_modis(
        modis_lai_ncdata_path,
        surface_space,
        start_date;
        time_interpolation_method = LinearInterpolation(),
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
    #    α_snow = Snow.ConstantAlbedoModel(FT(0.7))
    # Set β = 0 in order to regain model without density dependence
    α_snow = Snow.ZenithAngleAlbedoModel(
        FT(0.64),
        FT(0.06),
        FT(2);
        β = FT(0.4),
        x0 = FT(0.2),
    )
    horz_degree_res =
        sum(ClimaLand.Domains.average_horizontal_resolution_degrees(domain)) / 2 # mean of resolution in latitude and longitude, in degrees
    scf = Snow.WuWuSnowCoverFractionModel(
        FT(0.08),
        FT(1.77),
        FT(1.0),
        horz_degree_res,
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

    Y, p, cds = initialize(land)

    ic_path = ClimaLand.Artifacts.soil_ic_2008_50m_path(; context = context)
    evaluate!(p.snow.T, atmos.T, t0)
    ClimaLand.set_snow_initial_conditions!(
        Y,
        p,
        surface_space,
        ic_path,
        land.snow.parameters,
    )

    Y.soilco2.C .= FT(0.000412) # set to atmospheric co2, mol co2 per mol air
    Y.canopy.hydraulics.ϑ_l.:1 .= plant_ν
    evaluate!(Y.canopy.energy.T, atmos.T, t0)
    T_bounds = extrema(Y.canopy.energy.T)

    ClimaLand.set_soil_initial_conditions!(
        Y,
        ν,
        θ_r,
        subsurface_space,
        ic_path,
        land.soil,
        T_bounds,
    )
    set_initial_cache! = make_set_initial_cache(land)
    exp_tendency! = make_exp_tendency(land)
    imp_tendency! = ClimaLand.make_imp_tendency(land)
    jacobian! = ClimaLand.make_jacobian(land)
    set_initial_cache!(p, Y, t0)

    # set up jacobian info
    jac_kwargs = (;
        jac_prototype = ClimaLand.FieldMatrixWithSolver(Y),
        Wfact = jacobian!,
    )

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

    updateat = [promote(t0:(ITime(3600 * 3)):tf...)...]
    drivers = ClimaLand.get_drivers(land)
    updatefunc = ClimaLand.make_update_drivers(drivers)

    # ClimaDiagnostics
    # num_points is the resolution of the output diagnostics
    # These are currently chosen to get a 1:1 ration with the number of
    # simulation points, ~101x101x4x4
    nc_writer = ClimaDiagnostics.Writers.NetCDFWriter(
        subsurface_space,
        outdir;
        start_date,
        num_points = (570, 285, 15),
    )

    diags = ClimaLand.default_diagnostics(
        land,
        start_date;
        output_writer = nc_writer,
        output_vars = :short,
        average_period = :monthly,
    )

    diagnostic_handler =
        ClimaDiagnostics.DiagnosticsHandler(diags, Y, p, t0; dt = Δt)

    diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler)

    driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)

    walltime_info = WallTimeInfo()
    every1000steps(u, t, integrator) = mod(integrator.step, 1000) == 0
    report = let wt = walltime_info
        (integrator) -> report_walltime(wt, integrator)
    end
    report_cb = SciMLBase.DiscreteCallback(every1000steps, report)

    mask = ClimaLand.landsea_mask(domain)
    nancheck_freq = Dates.Month(1)
    nancheck_cb = ClimaLand.NaNCheckCallback(
        nancheck_freq,
        start_date,
        t0,
        Δt;
        mask = mask,
    )

    return prob,
    SciMLBase.CallbackSet(driver_cb, diag_cb, report_cb)
end

function setup_and_solve_problem(; greet = false)

    t0 = 0.0
    seconds = 1.0
    minutes = 60seconds
    hours = 60minutes # hours in seconds
    days = 24hours # days in seconds
    years = 366days # years in seconds - 366 to make sure we capture at least full years
    # 10 years in seconds for very long run and 2 years in seconds otherwise
    tf = LONGER_RUN ? 1years : 2years
    Δt = 450.0
    start_date = LONGER_RUN ? DateTime(2004) : DateTime(2008)
    nelements = (101, 15)
    if greet
        @info "Run: Global Soil-Canopy-Snow Model"
        @info "Resolution: $nelements"
        @info "Timestep: $Δt s"
        @info "Duration: $(tf - t0) s"
    end

    t0 = ITime(t0, epoch = start_date)
    tf = ITime(tf, epoch = start_date)
    Δt = ITime(Δt, epoch = start_date)
    t0, tf, Δt = promote(t0, tf, Δt)
    prob, cb = setup_prob(t0, tf, Δt, start_date; nelements)

    # Define timestepper and ODE algorithm
    stepper = CTS.ARS111()
    ode_algo = CTS.IMEXAlgorithm(
        stepper,
        CTS.NewtonsMethod(
            max_iters = 3,
            update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        ),
    )
    SciMLBase.solve(prob, ode_algo; dt = Δt, callback = cb, adaptive = false)
    return nothing
end

setup_and_solve_problem(; greet = true);
# read in diagnostics and make some plots!
#### ClimaAnalysis ####
simdir = ClimaAnalysis.SimDir(outdir)
short_names = [
    "gpp",
    "swc",
    "et",
    "shf",
    "swu",
    "lwu",
    "swe",
    "si",
    "lwp",
    "iwc",
    "trans",
    "msf",
    "snowc",
]

include(
    joinpath(
        pkgdir(ClimaLand),
        "experiments",
        "long_runs",
        "figures_function.jl",
    ),
)
make_figures(root_path, outdir, short_names)

include("leaderboard/leaderboard.jl")
diagnostics_folder_path = outdir
leaderboard_base_path = root_path
for data_source in ("ERA5", "ILAMB")
    compute_monthly_leaderboard(
        leaderboard_base_path,
        diagnostics_folder_path,
        data_source,
    )
    compute_seasonal_leaderboard(
        leaderboard_base_path,
        diagnostics_folder_path,
        data_source,
    )
end

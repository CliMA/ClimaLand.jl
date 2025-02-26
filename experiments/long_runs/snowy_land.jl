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
using ClimaCore
using ClimaUtilities.ClimaArtifacts
import ClimaUtilities.OnlineLogging: WallTimeInfo, report_walltime
import ClimaUtilities.TimeManager: ITime, date

import ClimaDiagnostics
import ClimaAnalysis
import ClimaAnalysis.Visualize as viz
import ClimaUtilities
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
time_interpolation_method =
    LONGER_RUN ? LinearInterpolation() : LinearInterpolation(PeriodicCalendar())
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
    forcing = ClimaLand.prescribed_forcing_era5(
        era5_ncdata_path,
        surface_space,
        start_date,
        earth_param_set,
        FT;
        time_interpolation_method = time_interpolation_method,
    )

    f_over = FT(3.28) # 1/m
    R_sb = FT(1.484e-4 / 1000) # m/s
    scalar_soil_params = (; f_over, R_sb)

    α_snow = FT(0.67)
    scalar_snow_params = (; α_snow,Δt)

    # Energy Balance model
    ac_canopy = FT(2.5e3)
    # Plant Hydraulics and general plant parameters
    K_sat_plant = FT(5e-9) # m/s # seems much too small?
    ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value is -4 MPa
    Weibull_param = FT(4) # unitless, Holtzman's original c param value
    a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
    plant_ν = FT(1.44e-4)
    plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
    h_leaf = FT(1.0) # canopy height in m
    scalar_canopy_params = (;
        ac_canopy,
        K_sat_plant,
        a,
        ψ63,
        Weibull_param,
        plant_ν,
        plant_S_s,
        h_leaf,
    )
    
    domain = ClimaLand.global_domain(FT; nelements = nelements, context = context)
    earth_param_set = LP.LandParameters(FT)
    land = global_land_model(FT,
                             scalar_soil_params,
                             scalar_canopy_params,
                             scalar_snow_params,
                             earth_param_set;
                             context = context,
                             domain = domain,
                             forcing = forcin
                             )
    
    Y, p, cds = initialize(land)

    ic_path = ClimaLand.Artifacts.soil_ic_2008_50m_path(; context = context)
    ClimaLand.set_soil_initial_conditions!(Y, ν, θ_r, subsurface_space, ic_path)
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
    # These are currently chosen to get a 1:1 ratio with the number of
    # simulation points, ~101x101x6
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

    mask = ClimaLand.landsea_mask(surface_space)
    nancheck_freq = Dates.Month(1)
    nancheck_cb = ClimaLand.NaNCheckCallback(
        nancheck_freq,
        start_date,
        t0,
        Δt;
        mask = mask,
    )

    return prob,
    SciMLBase.CallbackSet(driver_cb, diag_cb, report_cb, nancheck_cb)
end

function setup_and_solve_problem(; greet = false)

    t0 = 0.0
    seconds = 1.0
    minutes = 60seconds
    hours = 60minutes # hours in seconds
    days = 24hours # days in seconds
    years = 366days # years in seconds - 366 to make sure we capture at least full years
    # 10 years in seconds for very long run and 2 years in seconds otherwise
    tf = LONGER_RUN ? 10years : 2years
    Δt = 450.0
    start_date = DateTime(2008)
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
compute_leaderboard(leaderboard_base_path, diagnostics_folder_path)

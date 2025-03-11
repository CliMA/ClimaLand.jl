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
import ClimaUtilities.TimeManager: ITime

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

    f_over = FT(3.28) # 1/m
    R_sb = FT(1.484e-4 / 1000) # m/s
    other_soil_params = (; f_over, R_sb)

    α_snow = FT(0.67)
    other_snow_params = (; α_snow,)

    # Energy Balance model
    ac_canopy = FT(2.5e3)
    # Plant Hydraulics and general plant parameters
    SAI = FT(0.0) # m2/m2
    f_root_to_shoot = FT(3.5)
    RAI = FT(1.0)
    K_sat_plant = FT(5e-9) # m/s # seems much too small?
    ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value is -4 MPa
    Weibull_param = FT(4) # unitless, Holtzman's original c param value
    a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
    plant_ν = FT(1.44e-4)
    plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
    n_stem = 0
    n_leaf = 1
    h_stem = FT(0.0)
    h_leaf = FT(1.0)
    zmax = FT(0.0)
    other_canopy_params = (;
        ac_canopy,
        SAI,
        f_root_to_shoot,
        RAI,
        K_sat_plant,
        a,
        ψ63,
        Weibull_param,
        plant_ν,
        plant_S_s,
        n_stem,
        n_leaf,
        h_stem,
        h_leaf,
        zmax,
    )

    land, Y, p, cds = ClimaLand.land_model_setup(
        FT;
        earth_param_set,
        context,
        nelements,
        start_date,
        t0,
        Δt,
        other_soil_params,
        other_canopy_params,
        other_snow_params,
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
    # These are currently chosen to get a 1:1 ratio with the number of
    # simulation points, ~101x101x6
    subsurface_space = land.soil.domain.space.subsurface
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

    nancheck_freq = Dates.Month(1)
    nancheck_cb = ClimaLand.NaNCheckCallback(nancheck_freq, start_date, t0, Δt)

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
    tf = 2years # 2 years in seconds
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
short_names_bio = ["gpp", "ct", "lwp"]
short_names_water = ["swc", "si", "sr", "swe"]
short_names_other = ["swu", "lwu", "et"]
group_names = ["bio", "water", "other"]
months_id = [1, 4, 7, 10]
for (group_id, group) in
    enumerate([short_names_bio, short_names_water, short_names_other])
    fig =
        CairoMakie.Figure(size = (600 * length(months_id), 300 * length(group)))
    for (var_id, short_name) in enumerate(group)
        var = get(simdir; short_name)
        times = ClimaAnalysis.times(var)
        CairoMakie.Label(
            fig[var_id, 0],
            short_name,
            tellheight = false,
            tellwidth = false,
            fontsize = 20,
        )
        for (t_id, t) in pairs(times[months_id])
            layout = fig[var_id, t_id] = CairoMakie.GridLayout()
            kwargs = ClimaAnalysis.has_altitude(var) ? Dict(:z => 1) : Dict()
            ClimaAnalysis.Visualize.heatmap2D_on_globe!(
                layout,
                ClimaAnalysis.slice(var, time = t; kwargs...),
                mask = ClimaAnalysis.Visualize.oceanmask(),
                more_kwargs = Dict(
                    :mask => ClimaAnalysis.Utils.kwargs(color = :white),
                ),
            )
        end
    end
    months = Dates.monthname.(1:12) .|> x -> x[1:3]
    for (idx, m_id) in enumerate(months_id)
        CairoMakie.Label(
            fig[0, idx],
            months[m_id],
            tellwidth = false,
            fontsize = 20,
        )
    end
    group_name = group_names[group_id]

    CairoMakie.save(joinpath(root_path, "$(group_name).png"), fig)
end

short_names =
    ["gpp", "swc", "et", "ct", "swe", "si", "lwp", "iwc", "trans", "msf"]

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

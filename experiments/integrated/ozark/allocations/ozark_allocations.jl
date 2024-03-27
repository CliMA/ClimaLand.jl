import SciMLBase
import ClimaTimeSteppers as CTS
using ClimaCore
using Plots
using Statistics
using Dates

using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaParams

PROFILING = false
try
    import Profile, ProfileCanvas
    global PROFILING = true
    @info "ProfileCanvas found, running with profiler"
catch
end
global FT = Float64
climaland_dir = pkgdir(ClimaLand)
global site_ID = "US-MOz"
savedir = joinpath(climaland_dir, "experiments/integrated/ozark/allocations")
global earth_param_set = LP.LandParameters(FT)

# Create model and setup some timestepping parameters
include(
    joinpath(
        climaland_dir,
        "experiments/integrated/ozark/ozark_simulation_setup.jl",
    ),
)
# Setup initial conditions
include(
    joinpath(
        climaland_dir,
        "experiments/integrated/ozark/ozark_initial_conditions.jl",
    ),
)
tf = t0 + dt * 300
updateat = Array(t0:dt:tf)
updatefunc = ClimaLand.make_update_drivers(atmos, radiation)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction(
        T_exp! = exp_tendency!,
        dss! = ClimaLand.dss!,
        T_imp! = nothing,
    ),
    Y,
    (t0, tf),
    p,
)

sol = SciMLBase.solve(
    prob,
    ode_algo;
    dt = dt,
    callback = driver_cb,
    adaptive = false,
)
if PROFILING
    # Now that we compiled, solve again but collect profiling information
    include(
        joinpath(
            climaland_dir,
            "experiments/integrated/ozark/ozark_initial_conditions.jl",
        ),
    )
    updateat = Array(t0:dt:tf)
    updatefunc = ClimaLand.make_update_drivers(atmos, radiation)
    driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
    prob = SciMLBase.ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = exp_tendency!,
            dss! = ClimaLand.dss!,
            T_imp! = nothing,
        ),
        Y,
        (t0, tf),
        p,
    )
    Profile.@profile SciMLBase.solve(
        prob,
        ode_algo;
        dt = dt,
        callback = driver_cb,
    )
    results = Profile.fetch()
    flame_file = joinpath(savedir, "flame.html")
    ProfileCanvas.html_file(flame_file, results)
    @info "Save compute flame to flame_file"
    include(
        joinpath(
            climaland_dir,
            "experiments/integrated/ozark/ozark_initial_conditions.jl",
        ),
    )
    updateat = Array(t0:dt:tf)
    updatefunc = ClimaLand.make_update_drivers(atmos, radiation)
    driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
    prob = SciMLBase.ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = exp_tendency!,
            dss! = ClimaLand.dss!,
            T_imp! = nothing,
        ),
        Y,
        (t0, tf),
        p,
    )
    Profile.Allocs.@profile sample_rate = 1.0 SciMLBase.solve(
        prob,
        ode_algo;
        dt = dt,
        callback = driver_cb,
    )
    results = Profile.Allocs.fetch()
    profile = ProfileCanvas.view_allocs(results)
    alloc_flame_file = joinpath(savedir, "alloc_flame.html")
    ProfileCanvas.html_file(alloc_flame_file, profile)
    @info "Save allocation flame to alloc_flame_file"
end

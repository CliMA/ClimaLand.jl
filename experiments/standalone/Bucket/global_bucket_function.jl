# # Global bucket run using analytic albedo

# The code sets up and runs the bucket for 7 days using analytic albedo as
# a function of space, and analytic atmospheric and radiative forcings.
# This is the simplest global bucket run we test, and it is used to verify
# that the model can run on both CPU and GPU devices, with only minor
# computational differences between the results.

# Outputs:
# The final state of the simulation is saved to a CSV file so we can compare
# between CPU and GPU runs.
# Plots of the temporal evolution of water content, snow cover fraction,
# surface temperature, evaporation, and surface energy flux.

import SciMLBase
import ClimaComms
ClimaComms.@import_required_backends
using CairoMakie
using Dates
using DelimitedFiles
using Statistics

import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
import ClimaUtilities.OutputPathGenerator: generate_output_path

import ClimaTimeSteppers as CTS
import NCDatasets
using ClimaCore
using ClimaCore: Remapping, Geometry
import ClimaComms
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaParams
using ClimaLand.Bucket:
    BucketModel, BucketModelParameters, PrescribedBaregroundAlbedo
using ClimaLand.Domains: coordinates, Column
using ClimaLand:
    initialize,
    make_update_aux,
    make_exp_tendency,
    make_set_initial_cache,
    PrescribedAtmosphere,
    PrescribedRadiativeFluxes

using ClimaDiagnostics
using ClimaAnalysis
import ClimaAnalysis.Visualize as viz
using ClimaUtilities

"""
   compute_extrema(v)

Computes and returns the minimum value in `v` and
the maximum value in `v`, as a tuple, assuming that
`v` is a vector of arrays.
"""
function compute_extrema(v)
    maxes = [maximum(u) for u in v]
    mins = [minimum(u) for u in v]
    return (minimum(mins), maximum(maxes))
end

anim_plots = false
FT = Float64;
context = ClimaComms.context()
ClimaComms.init(context)
earth_param_set = LP.LandParameters(FT);
# Use separate output directory for CPU and GPU runs to avoid race condition
device_suffix =
    typeof(context.device) <: ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
outdir = "experiments/standalone/Bucket/artifacts_function_$(device_suffix)"

# Construct simulation domain
soil_depth = FT(3.5);
bucket_domain = ClimaLand.Domains.SphericalShell(;
    radius = FT(6.3781e6),
    depth = soil_depth,
    nelements = (50, 10),
    npolynomial = 1,
    dz_tuple = FT.((1.0, 0.05)),
);
surface_space = bucket_domain.space.surface

# Set up parameters
α_bareground_func = (coordinate_point) -> 0.2;
α_snow = FT(0.8);
albedo =
    PrescribedBaregroundAlbedo{FT}(α_snow, α_bareground_func, surface_space);
σS_c = FT(0.2);
W_f = FT(0.15);
z_0m = FT(1e-2);
z_0b = FT(1e-3);
κ_soil = FT(0.7);
ρc_soil = FT(2e6);
τc = FT(3600);
t0 = 0.0;
tf = 7 * 86400;
Δt = 3600.0;

bucket_parameters = BucketModelParameters(FT; albedo, z_0m, z_0b, τc);
start_date = DateTime(2005);

# Precipitation:
precip = (t) -> 0;
snow_precip = (t) -> -5e-7 * (t < 1 * 86400);
# Diurnal temperature variations:
T_atmos = (t) -> 275.0 + 5.0 * sin(2.0 * π * t / 86400 - π / 2);
# Constant otherwise:
u_atmos = (t) -> 3.0;
q_atmos = (t) -> 0.001;
h_atmos = FT(2);
P_atmos = (t) -> 101325;

bucket_atmos = PrescribedAtmosphere(
    TimeVaryingInput(precip),
    TimeVaryingInput(snow_precip),
    TimeVaryingInput(T_atmos),
    TimeVaryingInput(u_atmos),
    TimeVaryingInput(q_atmos),
    TimeVaryingInput(P_atmos),
    start_date,
    h_atmos,
    earth_param_set,
);

# Prescribed radiation -- a prescribed downwelling SW diurnal cycle, with a
# peak at local noon, and a prescribed downwelling LW radiative
# flux, assuming the air temperature is on average 275 degrees
# K with a diurnal amplitude of 5 degrees K:
SW_d = (t) -> max(1361 * sin(2π * t / 86400 - π / 2), 0.0);
LW_d = (t) -> 5.67e-8 * (275.0 + 5.0 * sin(2.0 * π * t / 86400 - π / 2))^4;
bucket_rad = PrescribedRadiativeFluxes(
    FT,
    TimeVaryingInput(SW_d),
    TimeVaryingInput(LW_d),
    start_date,
);

model = BucketModel(
    parameters = bucket_parameters,
    domain = bucket_domain,
    atmosphere = bucket_atmos,
    radiation = bucket_rad,
);

Y, p, coords = initialize(model);

Y.bucket.T .= FT(270);
Y.bucket.W .= FT(0.05);
Y.bucket.Ws .= FT(0.0);
Y.bucket.σS .= FT(0.08);

set_initial_cache! = make_set_initial_cache(model);
set_initial_cache!(p, Y, t0);
exp_tendency! = make_exp_tendency(model);
timestepper = CTS.RK4()
ode_algo = CTS.ExplicitAlgorithm(timestepper)
prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction(T_exp! = exp_tendency!, dss! = ClimaLand.dss!),
    Y,
    (t0, tf),
    p,
);

# ClimaDiagnostics
output_dir = ClimaUtilities.OutputPathGenerator.generate_output_path(outdir)

space = bucket_domain.space.subsurface

nc_writer = ClimaDiagnostics.Writers.NetCDFWriter(space, output_dir; start_date)

diags =
    ClimaLand.default_diagnostics(model, start_date; output_writer = nc_writer)

diagnostic_handler =
    ClimaDiagnostics.DiagnosticsHandler(diags, Y, p, t0; dt = Δt)

diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler)

updateat = collect(t0:(Δt * 3):tf);
drivers = ClimaLand.get_drivers(model)
updatefunc = ClimaLand.make_update_drivers(drivers)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)

sol = ClimaComms.@time ClimaComms.device() SciMLBase.solve(
    prob,
    ode_algo;
    dt = Δt,
    callback = SciMLBase.CallbackSet(driver_cb, diag_cb),
)

#### ClimaAnalysis ####

# all
simdir = ClimaAnalysis.SimDir(output_dir)
short_names_2D = [
    "swa",
    "rn",
    "tsfc",
    "qsfc",
    "lhf",
    "rae",
    "shf",
    "vflux",
    "rhosfc",
    "wsoil",
    "wsfc",
    "ssfc",
]
short_names_3D = ["tsoil"]
for short_name in vcat(short_names_2D..., short_names_3D...)
    var = get(simdir; short_name)
    fig = CairoMakie.Figure(size = (800, 600))
    kwargs = short_name in short_names_2D ? Dict() : Dict(:z => 1)
    viz.plot!(fig, var, lon = 0, lat = 0; kwargs...)
    CairoMakie.save(joinpath(output_dir, "$short_name.png"), fig)
end

# Interpolate to grid
space = axes(coords.surface)
longpts = range(-180.0, 180.0, 21)
latpts = range(-90.0, 90.0, 21)
hcoords = [Geometry.LatLongPoint(lat, long) for long in longpts, lat in latpts]
remapper = Remapping.Remapper(space, hcoords)

W = Array(Remapping.interpolate(remapper, sol.u[end].bucket.W))
Ws = Array(Remapping.interpolate(remapper, sol.u[end].bucket.Ws))
σS = Array(Remapping.interpolate(remapper, sol.u[end].bucket.σS))
T_sfc = Array(Remapping.interpolate(remapper, prob.p.bucket.T_sfc))


# save prognostic state to CSV - for comparison between
# GPU and CPU output
open(joinpath(output_dir, "tf_state_$(device_suffix)_function.txt"), "w") do io
    writedlm(io, hcat(T_sfc[:], W[:], Ws[:], σS[:]), ',')
end;

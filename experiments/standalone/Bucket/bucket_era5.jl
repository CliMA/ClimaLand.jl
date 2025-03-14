# # Global/regional bucket run using spatial map albedo

# The code sets up and runs the bucket for 7 days using albedo read in from
# a file containing static data over the globe or on a region, and ERA5 atmospheric and
# radiative forcings.
# Moving forward, this driver will serve as our more complex global bucket run,
# eventually running for a longer time period (1+ year) and using temporally
# varying atmospheric and radiative forcing data.
# This driver is used to verify that this more complex version of the model can
# run on both CPU and GPU, with only minor computational differences between the results.

# Outputs:
# The final state of the simulation is saved to a CSV file so we can compare
# between CPU and GPU runs.
# Plots of the temporal evolution of water content, snow cover fraction,
# surface temperature, evaporation, and surface energy flux.
import SciMLBase
import ClimaComms
ClimaComms.@import_required_backends
using Dates
using DelimitedFiles
using Statistics
using ClimaDiagnostics
using ClimaUtilities.ClimaArtifacts
import Interpolations
import ClimaAnalysis
import ClimaAnalysis
import GeoMakie
import ClimaAnalysis.Visualize as viz
using CairoMakie
import ClimaUtilities
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.OutputPathGenerator: generate_output_path
import ClimaUtilities.TimeManager: ITime
import ClimaTimeSteppers as CTS
import NCDatasets
using ClimaCore
using ClimaCore: Remapping, Geometry
import ClimaParams as CP
import ClimaComms
import ClimaLand
import ClimaLand.Parameters as LP
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

function compute_clims(v)
    means = [mean(u) for u in v]
    sigmas = [std(u) for u in v]
    return (minimum(means) - maximum(sigmas), maximum(means) + maximum(sigmas))
end

# Set to true if you want to run a regional simulation. By default, it is false,
# unless the `CLIMALAND_CI_REGIONAL_BUCKET` environment variable is defined.
regional_simulation = haskey(ENV, "CLIMALAND_CI_REGIONAL_BUCKET")
regional_str = regional_simulation ? "_regional" : ""
time_interpolation_method = LinearInterpolation(PeriodicCalendar())
regridder_type = :InterpolationsRegridder
FT = Float64;
context = ClimaComms.context()
ClimaComms.init(context)
earth_param_set = LP.LandParameters(FT);
# Use separate output directory for CPU and GPU runs to avoid race condition
device_suffix =
    typeof(context.device) <: ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
outdir = "experiments/standalone/Bucket/artifacts_era5$(regional_str)_$(device_suffix)"
output_dir = ClimaUtilities.OutputPathGenerator.generate_output_path(outdir)

# Set up simulation domain
soil_depth = FT(3.5);
if regional_simulation
    @info "Running regional simulation"
    center_long, center_lat = FT(-118.14452), FT(34.14778)
    # Half size of the grid in meters
    delta_m = FT(50_000)
    bucket_domain = ClimaLand.Domains.HybridBox(;
        xlim = (delta_m, delta_m),
        ylim = (delta_m, delta_m),
        zlim = (-soil_depth, FT(0.0)),
        longlat = (center_long, center_lat),
        nelements = (30, 30, 10),
        npolynomial = 1,
    )
else
    bucket_domain = ClimaLand.Domains.SphericalShell(;
        radius = FT(6.3781e6),
        depth = soil_depth,
        nelements = (30, 10),
        npolynomial = 1,
        dz_tuple = FT.((1.0, 0.05)),
    )
end
start_date = DateTime(2008);

# Set up parameters
σS_c = FT(0.2);
W_f = FT(0.15);
z_0m = FT(1e-2);
z_0b = FT(1e-3);
κ_soil = FT(0.7);
ρc_soil = FT(2e6);
τc = FT(3600);
t0 = 0.0;
tf = 14 * 86400;
Δt = 3600.0 / 3;

t0 = ITime(t0, epoch = start_date)
tf = ITime(tf, epoch = start_date)
Δt = ITime(Δt, epoch = start_date)
t0, tf, Δt = promote(t0, tf, Δt)

# Construct albedo parameter object using static map
surface_space = bucket_domain.space.surface
subsurface_space = bucket_domain.space.subsurface
α_snow = FT(0.8)
albedo = PrescribedBaregroundAlbedo{FT}(α_snow, surface_space);

bucket_parameters = BucketModelParameters(FT; albedo, z_0m, z_0b, τc);

# Forcing data
era5_ncdata_path =
    ClimaLand.Artifacts.era5_land_forcing_data2008_folder_path(; context)
bucket_atmos, bucket_rad = ClimaLand.prescribed_forcing_era5(
    era5_ncdata_path,
    surface_space,
    start_date,
    earth_param_set,
    FT;
    time_interpolation_method = time_interpolation_method,
    regridder_type = regridder_type,
)

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

saveat = [promote(t0:(3 * Δt):tf...)...];
saved_values = (;
    t = Array{typeof(t0)}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
);
saving_cb = ClimaLand.NonInterpSavingCallback(saved_values, saveat);
updateat = copy(saveat)
drivers = ClimaLand.get_drivers(model)
updatefunc = ClimaLand.make_update_drivers(drivers)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)

# Diagnostics
nc_writer = ClimaDiagnostics.Writers.NetCDFWriter(
    subsurface_space,
    output_dir;
    start_date,
)
diags = ClimaLand.default_diagnostics(
    model,
    start_date;
    output_writer = nc_writer,
    average_period = :daily,
)

diagnostic_handler =
    ClimaDiagnostics.DiagnosticsHandler(diags, Y, p, t0; dt = Δt)

diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler)
cb = SciMLBase.CallbackSet(driver_cb, saving_cb, diag_cb)

sol = ClimaComms.@time ClimaComms.device() SciMLBase.solve(
    prob,
    ode_algo;
    dt = Δt,
    saveat = saveat,
    callback = cb,
);

simdir = ClimaAnalysis.SimDir(output_dir)
short_names = ["rn", "tsfc", "qsfc", "lhf", "shf", "wsoil", "wsfc", "ssfc"]
for short_name in short_names
    var = get(simdir; short_name)
    t = ClimaAnalysis.times(var)[end]
    var = get(simdir; short_name)
    fig = CairoMakie.Figure(size = (800, 600))
    kwargs = ClimaAnalysis.has_altitude(var) ? Dict(:z => 1) : Dict()
    viz.heatmap2D_on_globe!(
        fig,
        ClimaAnalysis.slice(var, time = t; kwargs...),
        mask = viz.oceanmask(),
        more_kwargs = Dict(:mask => ClimaAnalysis.Utils.kwargs(color = :white)),
    )
    CairoMakie.save(joinpath(output_dir, "$(short_name)_$t.png"), fig)
end

# Interpolate to grid
space = axes(coords.surface)
num_pts = 21
if regional_simulation
    radius_earth = FT(6.378e6)
    xlim = (
        center_long - delta_m / (2radius_earth),
        center_long + delta_m / (2radius_earth),
    )
    ylim = (
        center_lat - delta_m / (2radius_earth),
        center_lat + delta_m / (2radius_earth),
    )
    longpts = range(xlim[1], xlim[2], num_pts)
    latpts = range(ylim[1], ylim[2], num_pts)

    # Temporary monkey patch until we use the new diagnostics.
    # This is used somewhere internally by the remapper.
    # TODO: Remove this and use ClimaDiagnostics instead
    ClimaCore.Geometry._coordinate(
        pt::ClimaCore.Geometry.LatLongPoint,
        ::Val{1},
    ) = ClimaCore.Geometry.LatPoint(pt.lat)
    ClimaCore.Geometry._coordinate(
        pt::ClimaCore.Geometry.LatLongPoint,
        ::Val{2},
    ) = ClimaCore.Geometry.LongPoint(pt.long)

else
    longpts = range(-180.0, 180.0, num_pts)
    latpts = range(-90.0, 90.0, num_pts)
    hcoords =
        [Geometry.LatLongPoint(lat, long) for long in longpts, lat in latpts]
end
hcoords = [Geometry.LatLongPoint(lat, long) for long in longpts, lat in latpts]
remapper = Remapping.Remapper(space, hcoords)

W = Array(Remapping.interpolate(remapper, sol.u[end].bucket.W))
Ws = Array(Remapping.interpolate(remapper, sol.u[end].bucket.Ws))
σS = Array(Remapping.interpolate(remapper, sol.u[end].bucket.σS))
T_sfc = Array(Remapping.interpolate(remapper, prob.p.bucket.T_sfc))

# save prognostic state to CSV - for comparison between GPU and CPU output
open(joinpath(output_dir, "tf_state_$(device_suffix)_era5.txt"), "w") do io
    writedlm(io, hcat(T_sfc[:], W[:], Ws[:], σS[:]), ',')
end;

# # Global bucket run using spatial map albedo

# The code sets up and runs the bucket for 7 days using albedo read in from
# a file containing static data over the globe, and analytic atmospheric and
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
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
using Dates
using DelimitedFiles
using Statistics
using ClimaUtilities.ClimaArtifacts
import Interpolations
import ClimaCoreMakie
using CairoMakie
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
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

anim_plots = true

# Set to true if you want to run a regional simulation. By default, it is false,
# unless the `CLIMALAND_CI_REGIONAL_BUCKET` environment variable is defined.
regional_simulation = haskey(ENV, "CLIMALAND_CI_REGIONAL_BUCKET")
regional_str = regional_simulation ? "_regional" : ""

regridder_type = :InterpolationsRegridder
FT = Float64;
context = ClimaComms.context()
earth_param_set = LP.LandParameters(FT);
outdir = joinpath(
    pkgdir(ClimaLand),
    "experiments/standalone/Bucket/artifacts_staticmap",
)
!ispath(outdir) && mkpath(outdir)

# Set up simulation domain
soil_depth = FT(3.5);
bucket_domain = ClimaLand.Domains.SphericalShell(;
    radius = FT(6.3781e6),
    depth = soil_depth,
    nelements = (30, 10),
    npolynomial = 1,
    dz_tuple = FT.((1.0, 0.05)),
);
ref_time = DateTime(2021);

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

# Construct albedo parameter object using static map
# Use separate regridding directory for CPU and GPU runs to avoid race condition
device_suffix =
    typeof(ClimaComms.context().device) <: ClimaComms.CPUSingleThreaded ?
    "cpu" : "gpu"
t_start = t0
surface_space = bucket_domain.space.surface
α_snow = FT(0.8)
albedo = PrescribedBaregroundAlbedo{FT}(α_snow, surface_space);

bucket_parameters = BucketModelParameters(FT; albedo, z_0m, z_0b, τc);

# Forcing data
era5_artifact_path =
    ClimaLand.Artifacts.era5_land_forcing_data2021_folder_path(; context)

# Below, the preprocess_func argument is used to
# 1. Convert precipitation to be negative (as it is downwards)
# 2. Convert accumulations over an hour to a rate per second
# Precipitation:
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
h_atmos = FT(10);


bucket_atmos = PrescribedAtmosphere(
    precip,
    snow_precip,
    T_atmos,
    u_atmos,
    q_atmos,
    P_atmos,
    ref_time,
    h_atmos,
    earth_param_set,
);

# Prescribed radiation -- a prescribed downwelling SW diurnal cycle, with a
# peak at local noon, and a prescribed downwelling LW radiative
# flux, assuming the air temperature is on average 275 degrees
# K with a diurnal amplitude of 5 degrees K:
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

bucket_rad = PrescribedRadiativeFluxes(FT, SW_d, LW_d, ref_time);

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

saveat = collect(t0:(Δt * 3):tf);
saved_values = (;
    t = Array{Float64}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
);
saving_cb = ClimaLand.NonInterpSavingCallback(saved_values, saveat);
updateat = copy(saveat)
drivers = ClimaLand.get_drivers(model)
updatefunc = ClimaLand.make_update_drivers(drivers)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
cb = SciMLBase.CallbackSet(driver_cb, saving_cb)

sol = ClimaComms.@time ClimaComms.device() SciMLBase.solve(
    prob,
    ode_algo;
    dt = Δt,
    saveat = saveat,
    callback = cb,
);

# Interpolate to grid
space = axes(coords.surface)
longpts = range(-180.0, 180.0, 21)
latpts = range(-90.0, 90.0, 21)
hcoords = [Geometry.LatLongPoint(lat, long) for long in longpts, lat in latpts]
remapper = Remapping.Remapper(space, hcoords)

W = [
    Array(Remapping.interpolate(remapper, sol.u[k].bucket.W)) for
    k in 1:length(sol.t)
];
Ws = [
    Array(Remapping.interpolate(remapper, sol.u[k].bucket.Ws)) for
    k in 1:length(sol.t)
];
σS = [
    Array(Remapping.interpolate(remapper, sol.u[k].bucket.σS)) for
    k in 1:length(sol.t)
];
T_sfc = [
    Array(
        Remapping.interpolate(remapper, saved_values.saveval[k].bucket.T_sfc),
    ) for k in 1:length(sol.t)
];
evaporation = [
    Array(
        Remapping.interpolate(
            remapper,
            saved_values.saveval[k].bucket.turbulent_fluxes.vapor_flux,
        ),
    ) for k in 1:length(sol.t)
];
F_sfc = [
    Array(
        Remapping.interpolate(
            remapper,
            saved_values.saveval[k].bucket.R_n .+
            saved_values.saveval[k].bucket.turbulent_fluxes.lhf .+
            saved_values.saveval[k].bucket.turbulent_fluxes.shf,
        ),
    ) for k in 1:length(sol.t)
];

sw_forcing = [
    Array(
        Remapping.interpolate(remapper, saved_values.saveval[k].drivers.SW_d),
    ) for k in 1:length(sol.t)
];

# save prognostic state to CSV - for comparison between GPU and CPU output
open(joinpath(outdir, "tf_state_$(device_suffix)_staticmap.txt"), "w") do io
    writedlm(io, hcat(T_sfc[end][:], W[end][:], Ws[end][:], σS[end][:]), ',')
end;
# animation settings
nframes = length(T_sfc) # hourly data
fig_ts = Figure(size = (1000, 1000))
for (i, (field_ts, field_name)) in enumerate(
    zip(
        [W, σS, T_sfc, evaporation, F_sfc, sw_forcing],
        ["W", "σS", "T_sfc", "evaporation", "F_sfc", "SW forcing"],
    ),
)
    if anim_plots
        fig = Figure(size = (1000, 1000))
        ax = Axis(
            fig[1, 1],
            xlabel = "Longitude",
            ylabel = "Latitude",
            title = field_name,
        )
        clims = compute_clims(field_ts)
        CairoMakie.Colorbar(fig[:, end + 1], colorrange = clims)
        outfile = joinpath(
            outdir,
            string("anim_$(device_suffix)_", field_name, ".mp4"),
        )
        record(
            fig,
            outfile,
            (nframes - 7 * 24):2:nframes;
            framerate = 3,
        ) do frame
            CairoMakie.heatmap!(
                longpts,
                latpts,
                field_ts[frame],
                colorrange = clims,
            )
        end
    end

end
outfile = joinpath(outdir, string("ts_$device_suffix.png"))
CairoMakie.save(outfile, fig_ts)

if device_suffix == "cpu"
    W_raw = sol.u[end].bucket.W
    σS_raw = sol.u[end].bucket.σS
    T_sfc_raw = saved_values.saveval[end].bucket.T_sfc
    fields = [W_raw, σS_raw, T_sfc_raw]
    titles = ["W", "σS", "T_sfc"]
    for (f, n) in zip(fields, titles)
        fig = Figure(size = (1000, 1000))
        ax = Axis(fig[1, 1])
        clims = extrema(f)
        ClimaCoreMakie.fieldheatmap!(ax, f, title = n)
        Colorbar(fig[:, end + 1], colorrange = clims)
        CairoMakie.save(
            joinpath(outdir, string(n, "raw_$(device_suffix).png")),
            fig,
        )
    end
end

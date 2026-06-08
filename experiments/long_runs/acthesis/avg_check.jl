
import ClimaComms
ClimaComms.@import_required_backends
using ClimaUtilities.ClimaArtifacts
import ClimaUtilities.TimeManager: ITime, date

import ClimaDiagnostics
import ClimaUtilities

import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP
using ClimaCore
using ClimaLand
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, solve!, step!

using CairoMakie, GeoMakie, ClimaAnalysis
import ClimaLand.LandSimVis as LandSimVis

using Dates, Flux, JLD2, StaticArrays, InteractiveUtils, Adapt, NCDatasets
NS = Base.get_extension(ClimaLand, :ConstrainedNeuralModelExt).NeuralSnow;

context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"

run_4748_path = "/home/acharbon/outputs/run_4748_outputs"
my_none_run = "/home/acharbon/outputs/none_global_monthly"
my_era5_rad = "/home/acharbon/outputs/era5_global/era_means_monthly_rad.nc"
clima_era5_rad = joinpath(
    ClimaLand.Artifacts.era5_monthly_averages_single_level_path(),
    "era5_monthly_averages_surface_single_level_197901-202410.nc",
)

my_era5_data = NCDataset(my_era5_rad)
clima_era5_data = NCDataset(clima_era5_rad)

function nanequal(x1, x2)
    if isnan(x1) && isnan(x2) #NaN == NaN usually returns false, need to override this
        return true
    else
        return x1 == x2
    end
end

function arr_nanequal(x1, x2)
    @assert size(x1) == size(x2)
    return count(nanequal.(x1, x2)) == prod(size(x1))
end

function nanmean(v)
    n = 0
    tot = 0.0
    for x in v
        if !isnan(x)
            n += 1
            tot += x
        end
    end
    n == 0 && return NaN
    return tot / n
end

function nanmean_weighted(v, w)
    n = 0.0
    tot = 0.0
    for (xi, wi) in zip(v, w)
        if !isnan(xi) && !isnan(wi)
            n += wi
            tot += xi * wi
        end
    end
    n == 0 && return NaN
    return tot / n
end

function nanmax(v)
    m = -Inf
    for x in v
        if !isnan(x) && x > m
            m = x
        end
    end
    return m
end

rmse(v1, v2) = sqrt(nanmean(abs2.(v1 .- v2)))

#Check my simulation outputs are identical to that of snowy_run
for v in ["lhf", "shf", "swu", "lwu", "swd", "lwd"]
    fname = "$(v)_1M_average.nc"
    clima_data = NCDataset(joinpath(run_4748_path, fname))
    my_data = NCDataset(joinpath(my_none_run, fname))
    @assert my_data[:date][1] == clima_data[:date][1]
    end_idx = length(clima_data[:date][:])
    clima_var = clima_data[v][1:end_idx, :, :]
    my_var = my_data[v][1:end_idx, :, :]
    @assert size(my_var) == size(clima_var)
    if arr_nanequal(my_var, clima_var)
        print("File $(v): MATCHES\n")
    else
        print("File $(v): DOES NOT MATCH\n")
    end
end

# Check our radiation files match at source level:
for v in [
    ("lhf", "mslhf"),
    ("shf", "msshf"),
    ("swu", "msuwswrf"),
    ("lwu", "msuwlwrf"),
]
    start_idx = findfirst(
        ==(Date(my_era5_data[:date][1] + Day(14))),
        Date.(clima_era5_data[:time][:]),
    )
    end_idx = findfirst(
        ==(Date(my_era5_data[:date][end] + Day(14))),
        Date.(clima_era5_data[:time][:]),
    )
    if v[1] == "swu"
        my_var = my_era5_data[:swd][:, :, :] .- my_era5_data[:swn][:, :, :]
    elseif v[1] == "lwu"
        my_var = my_era5_data[:lwd][:, :, :] .- my_era5_data[:lwn][:, :, :]
    else
        my_var = my_era5_data[v[1]][:, :, :]
    end
    clima_var = clima_era5_data[v[2]][:, :, start_idx:end_idx]
    clima_var_modified = circshift(clima_var[:, 1:(end - 1), :], (180, 0, 0)) #get rid of extra latitude index, shift longitude
    @assert size(my_var) == size(clima_var_modified)
    if arr_nanequal(my_var, clima_var_modified)
        print("File $(v): MATCHES\n")
    else
        print("File $(v): DOES NOT MATCH\n")
    end
end

# Check our radiation files match at the point at which they are actually compared to the simulation data
LVSE = Base.get_extension(ClimaLand, :LandSimulationVisualizationExt)
sim_var_dict = LVSE.get_sim_var_dict(my_none_run)
obs_var_dict = LVSE.get_obs_var_dict("ERA5")
mask_dict = LVSE.get_mask_dict("ERA5")
short_names = intersect(keys(sim_var_dict), keys(obs_var_dict))
sim_obs_comparsion_dict = Dict()

#Taken directly from leaderboard.jl:
for short_name in short_names
    @info short_name
    # Simulation data
    sim_var = sim_var_dict[short_name]()

    # Observational data
    obs_var = obs_var_dict[short_name](sim_var.attributes["start_date"])

    # Remove first spin_up_months months from simulation
    spin_up_months = 0 # keep it all for now
    spinup_cutoff = spin_up_months * 30 * 86400.0
    ClimaAnalysis.times(sim_var)[end] >= spinup_cutoff &&
        (sim_var = ClimaAnalysis.window(sim_var, "time", left = spinup_cutoff))

    # Get the first valid time and last valid time
    obs_times = ClimaAnalysis.times(obs_var)
    sim_times = ClimaAnalysis.times(sim_var)
    min_time = maximum(first.((obs_times, sim_times)))
    max_time = minimum(last.((obs_times, sim_times)))

    # Window OutputVars to restrain the times to those that are the same between
    # both OutputVars
    sim_var =
        ClimaAnalysis.window(sim_var, "time", left = min_time, right = max_time)
    obs_var =
        ClimaAnalysis.window(obs_var, "time", left = min_time, right = max_time)

    obs_var = ClimaAnalysis.shift_longitude(obs_var, -180.0, 180.0)
    obs_var = ClimaAnalysis.resampled_as(obs_var, sim_var)

    sim_obs_comparsion_dict[short_name] = (sim_var, obs_var)
end

for v in ["lhf", "shf", "swu", "lwu"]
    sim_var, obs_var = sim_obs_comparsion_dict[v]
    start_idx = findfirst(
        ==(DateTime(obs_var.attributes["start_date"]) + Hour(6)),
        my_era5_data[:date][:],
    )
    end_idx = findfirst(
        ==(
            Second(obs_var.dims["time"][end]) +
            DateTime(obs_var.attributes["start_date"]) +
            Hour(6),
        ),
        my_era5_data[:date][:],
    )
    clima_compare = obs_var.data
    if v == "swu"
        my_var = my_era5_data[:swd][:, :, :] .- my_era5_data[:swn][:, :, :]
    elseif v == "lwu"
        my_var = my_era5_data[:lwd][:, :, :] .- my_era5_data[:lwn][:, :, :]
    else
        my_var = my_era5_data[v][:, :, :]
    end
    my_compare = permutedims(my_var, (3, 1, 2))[start_idx:end_idx, :, :]
    if v in ["lhf", "shf"]
        my_compare = -1 * my_compare
    end
    print(
        "LARGEST DIFFERENCE BETWEEN $(v): ",
        nanmax(abs.(my_compare .- clima_compare)),
        "\n",
    )
end

#All are maximally off by 1e-12, meaning values compared are not the issue, its the aggregation that's the issue. Keep going:

#Check if masks are the same - for these sims the mask is just ClimaAnalysis.apply_oceanmask()
threshold = 0.5 #this is the default in ClimaAnalysis, even though we used 0.99 for the simulation:
cmask = ClimaAnalysis.generate_ocean_mask(NaN, 1.0, threshold = threshold)
temp_sim, temp_obs = sim_obs_comparsion_dict["lhf"]
m = Matrix{Float64}((!).(isnan.(cmask(temp_obs).data[1, :, :])))
my_mask = Matrix{Float64}(
    JLD2.load(
        joinpath("/home/acharbon/outputs/all_global_daily/landmask.jld2"),
    )["mask"],
)
print("MASK EQUIVALENCE: $(m == my_mask)")

#Check the RMSE / BIAS calc themselves

for short_name in ["lhf", "shf", "swu", "lwu"]
    print("EVALUATING $(short_name)\n")
    sim_var, obs_var = sim_obs_comparsion_dict[short_name]
    times = ClimaAnalysis.times(sim_var)
    mask = mask_dict[short_name](sim_var, obs_var)
    clima_rmse_vec = [
        ClimaAnalysis.global_rmse(
            ClimaAnalysis.slice(sim_var, time = t),
            ClimaAnalysis.slice(obs_var, time = t),
            mask = mask,
        ) for t in times
    ]
    clima_bias_vec = [
        ClimaAnalysis.global_bias(
            ClimaAnalysis.slice(sim_var, time = t),
            ClimaAnalysis.slice(obs_var, time = t),
            mask = mask,
        ) for t in times
    ]

    #rmse(sim, obs) = sqrt(nanmean(abs2.(sim .- obs)))
    #bias(sim, obs) = nanmean(sim .- obs)
    rmse_weighted(sim, obs, weights) =
        sqrt(nanmean_weighted(abs2.(sim .- obs), weights))
    bias_weighted(sim, obs, weights) = nanmean_weighted(sim .- obs, weights)

    my_mask[my_mask .== 0] .= NaN
    lat_weights = collect(transpose(hcat([cosd.((-90:89)) for _ in 1:360]...)))
    my_mask_shaped = reshape(my_mask, 1, size(my_mask)...)

    fname = "$(short_name)_1M_average.nc"
    my_sim = NCDataset(joinpath(my_none_run, fname))
    start_idx =
        findfirst(==(my_sim[:date][1] + Hour(6)), my_era5_data[:date][:])
    end_idx =
        findfirst(==(my_sim[:date][end] + Hour(6)), my_era5_data[:date][:])
    my_var = my_sim[short_name][:, :, :]
    if short_name == "swu"
        my_dat = my_era5_data[:swd][:, :, :] .- my_era5_data[:swn][:, :, :]
    elseif short_name == "lwu"
        my_dat = my_era5_data[:lwd][:, :, :] .- my_era5_data[:lwn][:, :, :]
    else
        my_dat = my_era5_data[short_name][:, :, :]
    end
    my_obs = permutedims(
        my_dat[:, :, start_idx:end_idx],
        (3, 1, 2),
    )
    if short_name in ["lhf", "shf"]
        my_obs = -1 * my_obs
    end
    my_var_masked = my_var .* my_mask_shaped
    my_obs_masked = my_obs .* my_mask_shaped
    n = size(my_var_masked)[1]

    my_rmse_vec = [
        rmse_weighted(
            my_var_masked[i, :, :],
            my_obs_masked[i, :, :],
            lat_weights,
        ) for i in 1:n
    ]
    my_bias_vec = [
        bias_weighted(
            my_var_masked[i, :, :],
            my_obs_masked[i, :, :],
            lat_weights,
        ) for i in 1:n
    ]

    rmse_diff = maximum(
        abs(t[1] - t[2]) for t in collect(zip(my_rmse_vec, clima_rmse_vec))
    )
    bias_diff = maximum(
        abs(t[1] - t[2]) for t in collect(zip(my_bias_vec, clima_bias_vec))
    )
    print(
        "$(short_name): MAXIMUM ABS DIFF BETWEEN RMSES: $(rmse_diff), (mean diff: $(nanmean(abs.(my_rmse_vec .- clima_rmse_vec))), mean size: $(nanmean(clima_rmse_vec)))\n",
    )
    print(
        "$(short_name): MAXIMUM ABS DIFF BETWEEN BIASES: $(bias_diff), (mean diff: $(nanmean(abs.(my_bias_vec .- clima_bias_vec))), mean size: $(nanmean(clima_bias_vec)))\n",
    )
end

#DIFFERENCES SO FAR:
# - you take out the first 12 months, I take out the first 3 months
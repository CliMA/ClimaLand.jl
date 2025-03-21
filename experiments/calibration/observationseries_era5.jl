using ClimaAnalysis
using ClimaLand
using ClimaLand.Artifacts
using ClimaCore
using Statistics
import EnsembleKalmanProcesses as EKP
using LinearAlgebra
using Dates
FT = Float64

group_by_years(years, data) = Dict(year => data[findall(y -> y == year, years)] for year in unique(years))

# Some variables different conventions/units between ERA5 and CliMA. These
# preprocessing_function transform everything into the CliMA conventions.
preprocessing_function = Dict("shf" => (-), "lhf" => (-))

# TODO: This won't be needed when ClimaLand shortvar follow CMIP naming convention
clima_to_era5 = Dict(
                     "lhf" => "mslhf",
                     "shf" => "msshf",
                     "swu" => "msuwswrf",
                     "lwu" => "msuwlwrf",
                    )

function year_from_seasonal_outputvar(out::ClimaAnalysis.OutputVar)
    # We add one Month because we want to map DJF to the year of JF
    years = year.(Month(1) .+ dates(out))
    unique!(years)
    if length(years) == 1
        return first(years)
    else
        error("Year is not well defined")
    end
end

# Load ERA5 data
era5_path = joinpath(
                     ClimaLand.Artifacts.era5_surface_data_1979_2024_path(),
                     "era5_monthly_averages_surface_single_level_197901-202410.nc",
                    )

# Get era5 data. We shift everything to be defined at the beginning of the month
# so that we can more easily ensure consistency with respect to CliMA data.
vars_temp = Dict(
                 var_name => OutputVar(era5_path, clima_to_era5[var_name], shift_by = Dates.firstdayofmonth) for var_name in variable_list
                )

# Get ClimaLand data (for regridding)
simdir = SimDir(
    "calibration_output_utki_sample/iteration_000/member_001/global_diagnostics/output_active/",
)

out_vars = Dict(var_name => get(simdir, var_name) for var_name in variable_list)

#= For developping
vars_temp = Dict(
var_name => window(var, "time", left=DateTime(2008), right=DateTime(2012)) for (var_name, var) in vars_temp
)
=#

# Resample variables to diagnostic grid
# note: this is very slow with 45 years of data (era5).
vars = Dict()
for var in keys(vars_temp)
    preprocessing_fun = get(preprocessing_function, var, identity)
    vars[var] = preprocessing_fun(ClimaAnalysis.resampled_as(
        vars_temp[var],
        out_vars[string(var)],
        dim_names = ("longitude", "latitude"),
    ))
end

# Extract seasonal data in a Dict named results, which contains:
# results[:all_seasons] is a dictionary that maps variables names into vector of size length(training_locations) of maps of years to vectors of size 4 - the seasonal average of the data. For example, results[:all_seasons]["lhf"][1][2016] is a vector of 4 elements, with the four seasonal averages for year 2016, for the first location.
# results[:seasonal_vars] is a dictionary that maps variables names into vectors of size length(training_locations), with elements that contain the variance across years for each season (length=4). For example, results[:seasonal_vars]["lhf"]
# results[:seasonal_means] is the same format as seasonal_vars, but contains mean instead of variance. For example, results[:seasonal_means]["lhf"]

results = Dict(k=>Dict(v => [] for v in variable_list) for k in (:all_seasons, :seasonal_vars, :seasonal_means))

"""
    fun_across_season(var_name, fun, lon, lat; split_fun = split_by_season)


Slice `var_name` across `lon` and `lat`, apply `split_fun`, extract the data as
`array`s, and apply `fun` over the resulting arrays.
"""
fun_across_season(var_name, fun, lon, lat; split_fun = split_by_season) = fun.(getproperty.(split_fun(slice(vars[var_name]; lon, lat)), :data))

years = year_from_seasonal_outputvar.(split_by_season_across_time(first(values(vars))))

years_to_delete = [y for y in unique(years) if y < year(start_date + spinup_period)]
valid_years = [y for y in unique(years) if !(y in years_to_delete)]

for (lon, lat) in training_locations
    seasonal = Dict(
                    var_name => group_by_years(years, fun_across_season(var_name, mean, lon, lat; split_fun = split_by_season_across_time)) for
                    var_name in variable_list
                   )

    # remove years before start_date + spinup_period
    foreach(var_name ->
            foreach(y->delete!(seasonal[var_name], y), years_to_delete),
             variable_list)

    variances = Dict(
                 var_name => fun_across_season(var_name, Statistics.var, lon, lat) for
        var_name in variable_list)

    variances = Dict(
                 var_name => fun_across_season(var_name, Statistics.var, lon, lat) for
        var_name in variable_list
    )

    means = Dict(
                 var_name => fun_across_season(var_name, mean, lon, lat) for
        var_name in variable_list
    )

    for var in variable_list
        push!(results[:all_seasons][var], seasonal[var])
        push!(results[:seasonal_vars][var], variances[var])
        push!(results[:seasonal_means][var], means[var])
    end
end

# Build noise
# Currently, replaces values in results[:seasonal_vars] with 25.
# For example, results[:seasonal_vars]["lhf"][1] is a 4 element vector of 25.0 values. (no annual variation)
# TODO: This should create a new Dict, maybe called noise
# This noise can be a combination of results[:seasonal_vars], results[:seasonal_means], weighted by latitudes, etc.
# For now, we just make it flat values (for example, for lhf, we assume the model error is 25.0 W m^-2 everywhere)

#noise = Dict()
#for var_name in variable_list
#    noise[var_name[i]] = Dict()
#end

# latitudes = map(x -> x[2], training_locations)
for var in values(results[:seasonal_vars])
    for i in 1:length(var)
        var[i] .= 5^2 # flat noise
        # IF scale by lat, # var[i] .= var[i] ./ max(cosd(lat), 0.1) # need to get lat
        # IF add model error, add results[:seasonal_means][var] .* 0.05 # or some factor
    end
end

# Generate observation series
obs_y = [
    EKP.combine_observations([
        EKP.combine_observations([
            EKP.Observation(
                Dict(
                    "samples" =>
                        results[:all_seasons][var_name][i][y],
                        "covariances" => Diagonal(results[:seasonal_vars][var_name][i]),
                    "names" => "$(var)_$(lon)_$(lat)_$(y)",
                ),
            ) for var_name in variable_list
        ]) for (i, (lon, lat)) in enumerate(training_locations)
    ]) for y in valid_years
]

n_samples = length(valid_years)
m_size = 1  # 1 year
given_batches =
    [collect(((i - 1) * m_size + 1):(i * m_size)) for i in 1:n_samples]
minibatcher = EKP.FixedMinibatcher(given_batches)
o_names = [string(y) for y in valid_years]

observationseries = EKP.ObservationSeries(obs_y, minibatcher, o_names)

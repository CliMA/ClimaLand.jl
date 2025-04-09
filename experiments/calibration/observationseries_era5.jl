using ClimaAnalysis
using ClimaLand
using ClimaLand.Artifacts
using ClimaCore
using Statistics
import EnsembleKalmanProcesses as EKP
using LinearAlgebra
using Dates
FT = Float64

var_list = [:lhf, :shf, :swu, :lwu]

preprocessing_function = Dict(:shf => (-), :lhf => (-))

# Move this function to ClimaAnalysis
function Dates.year(out::ClimaAnalysis.OutputVar)
    dates = Month(1) .+ Second.(times(out)) .+ DateTime(out.attributes["start_date"])
    years = year.(dates)
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

# Get era5 data
vars_temp = Dict(
    var => OutputVar(era5_path, varname, shift_by = Dates.firstdayofmonth) for (var, varname) in zip(
        var_list,
        ["mslhf", "msshf", "msuwswrf", "msuwlwrf"],
    )
)

# Get ClimaLand data (for regridding)
simdir = SimDir(
    "calibration_output_utki_sample/iteration_000/member_001/global_diagnostics/output_active/",
)
out_vars = Dict(var => get(simdir, var) for var in ["lhf", "swu", "lwu", "shf"])

# Resample variables to diagnostic grid
# note: this is very slow. is it normal?
vars = Dict()
for var in keys(vars_temp)
    vars[var] = ClimaAnalysis.resampled_as(
        vars_temp[var],
        out_vars[string(var)],
        dim_names = "longitude",
    )
end

# Extract seasonal data
results = Dict(k=>Dict(v => [] for v in var_list) for k in (:all_seasons, :seasonal_vars, :seasonal_means))

for (lon, lat) in training_locations

    seasonal_slices = Dict(
                           var => split_by_season_across_time(slice(vars[var], longitude = lon, latitude = lat)) for
                           var in var_list
                          )

#    seasonal_years = Dict(
#                          var => seasonal_slices[var][i].data

    seasonal = Dict(
                    var => [mean(seasonal_slices[var][i].data) for i in 1:length(seasonal_slices[var])]
                    for var in keys(seasonal_slices)
                   ) # Note: the edge season (1979 or 2024) might average over less than 3 months

    for var in var_list
        maybe_preprocess = get(preprocessing_function, var, identity)
        seasonal[var] .= maybe_preprocess.(seasonal[var])
    end

# OLD code below, with manual indexing instead of ClimaAnalysis function split_by_season_across_time
#    slices = Dict(
#        var => slice(vars[var], longitude = lon, latitude = lat).data for
#        var in [:lhf, :shf, :swu, :lwu]
#    )
#    slices[:lhf] .*= -1
#    slices[:shf] .*= -1
#    seasonal = Dict(
#        var => [mean(data[i:(i + 2)]) for i in 1:3:(length(data) - 2)] for
#        (var, data) in slices
#    )
    variances = Dict(
        v => [FT(var(seasonal[v][i:4:end])) for i in 1:4] for
        v in keys(seasonal)
    )
    means = Dict(
        var => [FT(mean(seasonal[var][i:4:end])) for i in 1:4] for
        var in keys(seasonal)
    )

    for var in var_list
        push!(results[:all_seasons][var], FT.(seasonal[var]))
        push!(results[:seasonal_vars][var], variances[var])
        push!(results[:seasonal_means][var], means[var])
    end
end

# latitudes = map(x -> x[2], training_locations)
for var in values(results[:seasonal_vars])
    for i in 1:length(var)
        var[i] .= 5^2 # flat noise
        # IF scale by lat, # var[i] .= var[i] ./ max(cosd(lat), 0.1) # need to get lat
        # IF add model error, add results[:seasonal_means][var] .* 0.05 # or some factor
    end
end

# Generate observation series
n_samples = 10
obs_y = [
    EKP.combine_observations([
        EKP.combine_observations([
            EKP.Observation(
                Dict(
                    "samples" =>
                        all_seasons[i][(1 + (4 * (y - 1))):(4 + (4 * (y - 1)))],
                    "covariances" => Diagonal(seasonal_vars[i]),
                    "names" => "$(var)_$(lon)_$(lat)_$(y)",
                ),
            ) for (var, all_seasons, seasonal_vars) in zip(
                ["lhf", "shf", "swu", "lwu"],
                [
                    results[:all_seasons][:lhf],
                    results[:all_seasons][:shf],
                    results[:all_seasons][:swu],
                    results[:all_seasons][:lwu],
                ],
                [
                    results[:seasonal_vars][:lhf],
                    results[:seasonal_vars][:shf],
                    results[:seasonal_vars][:swu],
                    results[:seasonal_vars][:lwu],
                ],
            )
        ]) for (i, (lon, lat)) in enumerate(training_locations)
    ]) for y in 31:(31 + n_samples)
]

m_size = 1  # 1 year
given_batches =
    [collect(((i - 1) * m_size + 1):(i * m_size)) for i in 1:n_samples]
minibatcher = EKP.FixedMinibatcher(given_batches)
o_names = [string(i) for i in 2009:(2009 + n_samples - 1)]

observationseries = EKP.ObservationSeries(obs_y, minibatcher, o_names)

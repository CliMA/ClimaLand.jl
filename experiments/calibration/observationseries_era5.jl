using ClimaAnalysis
using ClimaLand
using ClimaLand.Artifacts
using ClimaCore
using Statistics
import EnsembleKalmanProcesses as EKP
using LinearAlgebra
FT = Float64

# Load ERA5 data
era5_path = joinpath(
    ClimaLand.Artifacts.era5_surface_data_1979_2024_path(),
    "era5_monthly_averages_surface_single_level_197901-202410.nc",
)

# Get era5 data
vars_temp = Dict(
    var => OutputVar(era5_path, varname) for (var, varname) in zip(
        [:lhf, :shf, :swu, :lwu],
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
results = Dict(
    :all_seasons => Dict(:lhf => [], :shf => [], :swu => [], :lwu => []),
    :seasonal_vars => Dict(:lhf => [], :shf => [], :swu => [], :lwu => []),
    :seasonal_means => Dict(:lhf => [], :shf => [], :swu => [], :lwu => []),
)

for (lon, lat) in training_locations
    slices = Dict(
        var => slice(vars[var], longitude = lon, latitude = lat).data for
        var in [:lhf, :shf, :swu, :lwu]
    )
    slices[:lhf] .*= -1
    slices[:shf] .*= -1
    seasonal = Dict(
        var => [mean(data[i:(i + 2)]) for i in 1:3:(length(data) - 2)] for
        (var, data) in slices
    )
    variances = Dict(
        v => [FT(var(seasonal[v][i:4:end])) for i in 1:4] for
        v in keys(seasonal)
    )
    means = Dict(
        var => [FT(mean(seasonal[var][i:4:end])) for i in 1:4] for
        var in keys(seasonal)
    )

    for var in [:lhf, :shf, :swu, :lwu]
        push!(results[:all_seasons][var], FT.(seasonal[var]))
        push!(results[:seasonal_vars][var], variances[var])
        push!(results[:seasonal_means][var], means[var])
    end
end

# Apply scale factor
latitudes = map(x -> x[2], training_locations)
scale_factor = [5^2 / max(cosd(lat), 0.1) for lat in latitudes]
for var in values(results[:seasonal_vars])
    var .= scale_factor
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
                    "covariances" => Diagonal([seasonal_vars[i]]),
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

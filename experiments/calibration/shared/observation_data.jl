import ClimaAnalysis

filename = joinpath(
        ClimaLand.Artifacts.era5_surface_data2008_path(),
        "era5_monthly_surface_fluxes_200801-200812.nc",
    )

era5_lhf = ClimaAnalysis.OutputVar(filename, "mslhf")
era5_shf = ClimaAnalysis.OutputVar(filename, "msshf")
era5_lwu = ClimaAnalysis.OutputVar(filename, "msuwlwrf")
era5_swu = ClimaAnalysis.OutputVar(filename, "msuwswrf")

observations_vecs = []
for var in [era5_lhf, era5_shf, era5_lwu, era5_swu]
    var_global_average =
    ClimaAnalysis.average_lon(
        ClimaAnalysis.average_lat( # should be weighted_average_lat
            ClimaAnalysis.apply_oceanmask(
                var,
            )
        )
    ).data
    push!(observations_vecs, var_global_average)
end
observations = vcat(observations_vecs...)
observations[1:24] .*= -1 # different sign convention for lhf and shf



# Read in the era5 datafile
#era5_ds = Dataset(
#    joinpath(
#        ClimaLand.Artifacts.era5_surface_data2008_path(),
#        "era5_monthly_surface_fluxes_200801-200812.nc",
#    ),
#)

# Make the ERA5 target
#ERA5_target = []
#close = (x, y) -> abs(x - y) < 5e-1
#for (lon, lat) in training_locations
#    # Fetch slices of lhf and shf era5 data from the era5 dataset
#    lat_ind, lon_ind = findall((x) -> close(x, lat), era5_ds["latitude"][:])[1],
#    findall((x) -> close(x, lon + 180), era5_ds["longitude"][:])[1]
#    lhf_loc = vec(era5_ds["mslhf"][lon_ind, lat_ind, :][:, end, :])
#    shf_loc = vec(era5_ds["msshf"][lon_ind, lat_ind, :][:, end, :])
#    lwu_loc = vec(era5_ds["msuwlwrf"][lon_ind, lat_ind, :][:, end, :])
#    swu_loc = vec(era5_ds["msuwswrf"][lon_ind, lat_ind, :][:, end, :])
#
#    # Create Observation objects for lhf and shf
#    lhf_ERA5 = EKP.Observation(
#        Dict(
#            "samples" => lhf_loc,
#            "covariances" => cov(lhf_loc) * EKP.I,
#            "names" => "lhf_$(lon)_$(lat)",
#        ),
#    )
#
#    shf_ERA5 = EKP.Observation(
#        Dict(
#            "samples" => shf_loc,
#            "covariances" => cov(shf_loc) * EKP.I,
#            "names" => "shf_$(lon)_$(lat)",
#        ),
#    )
#
#    lwu_ERA5 = EKP.Observation(
#        Dict(
#            "samples" => lwu_loc,
#            "covariances" => cov(lwu_loc) * EKP.I,
#            "names" => "lwu_$(lon)_$(lat)",
#        ),
#    )
#    swu_ERA5 = EKP.Observation(
#        Dict(
#            "samples" => swu_loc,
#            "covariances" => cov(swu_loc) * EKP.I,
#            "names" => "swu_$(lon)_$(lat)",
#        ),
#    )
#
#    # Add the observations to the target
#    push!(ERA5_target, lhf_ERA5)
#    push!(ERA5_target, shf_ERA5)
#    push!(ERA5_target, lwu_ERA5)
#    push!(ERA5_target, swu_ERA5)
#end
#
#observations = EKP.get_obs(EKP.combine_observations(ERA5_target))

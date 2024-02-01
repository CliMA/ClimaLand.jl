import ArtifactWrappers as AW

export era5_t2m_sp_u10n_dataset_path, era5_t2m_sp_u10n_static_dataset_path

"""
    era5_t2m_sp_u10n_dataset_path()

Triggers the download of a subset of the global ERA5 dataset, if not
already downloaded, using Julia Artifacts, and returns the path to
this file.

This dataset includes temperature at 2m (`t2m`), surface pressure (`sp`), and
the u component of neutral wind at 10m (`u10n`). This data starts on 1/1/2021,
is sampled hourly, and covers 24 hours.
"""
function era5_t2m_sp_u10n_dataset_path()
    era5_t2m_sp_u10n_dataset = AW.ArtifactWrapper(
        @__DIR__,
        "era5_t2m_sp_u10n",
        AW.ArtifactFile[AW.ArtifactFile(
            url = "https://caltech.box.com/shared/static/zp8pv4c7y1lhch9s5hfr9ffi86ifsurl.nc",
            filename = "era5_t2m_sp_u10n_20210101.nc",
        ),],
    )
    path = AW.get_data_folder(era5_t2m_sp_u10n_dataset)
    return joinpath(path, "era5_t2m_sp_u10n_20210101.nc")
end

"""
    era5_t2m_sp_u10n_static_dataset_path()

Triggers the download of a subset of the global ERA5 dataset, if not
already downloaded, using Julia Artifacts, and returns the path to
this file.

This dataset includes temperature at 2m (`t2m`), surface pressure (`sp`), and
the u component of neutral wind at 10m (`u10n`). This data contains data on
1/1/2021, and does not have a time dimension.
"""
function era5_t2m_sp_u10n_static_dataset_path()
    era5_t2m_sp_u10n_static_dataset = AW.ArtifactWrapper(
        @__DIR__,
        "era5_t2m_sp_u10n_static",
        AW.ArtifactFile[AW.ArtifactFile(
            url = "https://caltech.box.com/shared/static/b9llj1mbe5ndktjyokl54kv8e7jslqhc.nc",
            filename = "era5_t2m_sp_u10n_20210101_static.nc",
        ),],
    )
    path = AW.get_data_folder(era5_t2m_sp_u10n_static_dataset)
    return joinpath(path, "era5_t2m_sp_u10n_20210101_static.nc")
end

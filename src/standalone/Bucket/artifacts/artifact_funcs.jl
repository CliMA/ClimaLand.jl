import ArtifactWrappers as AW

export cesm2_albedo_dataset_path, bareground_albedo_dataset_path

"""
    cesm2_albedo_dataset_path()

Triggers the download of the CESM2 land albedo dataset, if not
already downloaded, using Julia Artifacts, and returns the path to
this file.

This dataset contains monthly albedo data from 15/01/1850
to 15/12/2014.
"""
function cesm2_albedo_dataset_path()
    land_albedo_dataset = AW.ArtifactWrapper(
        @__DIR__,
        "land_albedo",
        AW.ArtifactFile[AW.ArtifactFile(
            url = "https://caltech.box.com/shared/static/h3fcq4hlfe0o4kaj7ghftv76tg7v879v.nc",
            filename = "esw_albedo_Amon_CESM2_historical_r1i1p1f1_gn_185001-201412_v2.nc",
        ),],
    )
    path = AW.get_data_folder(land_albedo_dataset)
    return joinpath(
        path,
        "esw_albedo_Amon_CESM2_historical_r1i1p1f1_gn_185001-201412_v2.nc",
    )
end

"""
    bareground_albedo_dataset_path()

Triggers the download of the average bareground land albedo dataset, if not
already downloaded, using Julia Artifacts, and returns the path to
this file.

This dataset does not contain a time component.
"""
function bareground_albedo_dataset_path()
    bareground_albedo_dataset = AW.ArtifactWrapper(
        @__DIR__,
        "bareground_albedo",
        AW.ArtifactFile[AW.ArtifactFile(
            url = "https://caltech.box.com/shared/static/ga9385kyl82t955dylsbnn4x51b412md.nc",
            filename = "bareground_albedo.nc",
        ),],
    )
    path = AW.get_data_folder(bareground_albedo_dataset)
    return joinpath(path, "bareground_albedo.nc")
end

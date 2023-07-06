import ArtifactWrappers as AW

export cesm2_albedo_dataset_path, bareground_albedo_dataset_path, mask_dataset_path

"""
    cesm2_albedo_dataset_path()

Triggers the download of the CESM2 land albedo dataset, if not
already downloaded, using Julia Artifacts, and returns the path to
this file.

This file contains surface albedo over the entire globe, and does include
some NaN values. It was created by Renato Braghiere by calculating
outgoing SW radiation / incoming SW radiation at the surface using
CESM2 data.
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

# https://caltech.box.com/shared/static/zx07ipsdkv8cbx9wmk4v24e49jkcyefq.nc
"""
    cesm2_land_albedo_dataset_path()

Triggers the download of the CESM2 LS3MIP land albedo dataset, if not
already downloaded, using Julia Artifacts, and returns the path to
this file.

This file contains surface albedo over land. NaN values have been replaced
by 1, since we expect that regions having 0 incoming radiation are
experiencing winter storms and therefore likely have high albedo at the
surface (snow, ice, etc). It was created by Renato Braghiere by
calculating outgoing SW radiation / incoming SW radiation at the surface
using CESM2 data.
"""
function cesm2_land_albedo_dataset_path()
    land_albedo_dataset = AW.ArtifactWrapper(
        @__DIR__,
        "land_albedo", # TODO change varname?
        AW.ArtifactFile[AW.ArtifactFile(
            url = "https://caltech.box.com/shared/static/zx07ipsdkv8cbx9wmk4v24e49jkcyefq.nc",
            filename = "rsalb_Amon_CESM2_land-hist_r1i1p1f1_gn_185001-201512.nc",
        ),],
    )
    path = AW.get_data_folder(land_albedo_dataset)
    return joinpath(
        path,
        "rsalb_Amon_CESM2_land-hist_r1i1p1f1_gn_185001-201512.nc",
    )
end
"""
    bareground_albedo_dataset_path()

Triggers the download of the average bareground land albedo dataset, if not
already downloaded, using Julia Artifacts, and returns the path to
this file.
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

"""
    mask_dataset_path()

Triggers the download of the land mask dataset, if not already downloaded,
and returns the path to this file.
"""
function mask_dataset_path()
    mask_dataset = AW.ArtifactWrapper(
        @__DIR__,
        "land_mask",
        AW.ArtifactFile[AW.ArtifactFile(
            url = "https://caltech.box.com/shared/static/vubmq84nhvbgdqayezguf3i1w6nqtwvu.ncc",
            filename = "seamask.nc",
        ),],
    )
    path = AW.get_data_folder(mask_dataset)
    return joinpath(path, "seamask.nc")
end

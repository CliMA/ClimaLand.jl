# All data is currently downloaded from Caltech box

### TUTORIAL DATA
# Data from Mizoguchi, M., 1990; used for Soil freezing front tutorial
# [Full citation:] Mizoguchi, M., 1990. Water, heat and salt transport in freezing soil.
#  Ph.D. thesis. (In Japanese.) University of Tokyo.
function mizoguchi_dataset_path()
    mizoguchi_dataset = AW.ArtifactWrapper(
        @__DIR__,
        "mizoguchi",
        AW.ArtifactFile[AW.ArtifactFile(
            url = "https://caltech.box.com/shared/static/3xbo4rlam8u390vmucc498cao6wmqlnd.csv",
            filename = "mizoguchi_all_data.csv",
        ),],
    )
    return AW.get_data_folder(mizoguchi_dataset)
end


### EXPERIMENT DATA
# Data from AmeriFlux Ozark FLUXNET site; used for Ozark FLUXNET experiment
# [Full citation:] Jeffrey Wood, Lianhong Gu (2021), AmeriFlux FLUXNET-1F US-MOz Missouri Ozark Site,
#  Ver. 3-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1854370
function ozark_fluxnet_dataset_path()
    ozark_fluxnet_dataset = AW.ArtifactWrapper(
        @__DIR__,
        "ameriflux_data-US-MOz",
        AW.ArtifactFile[AW.ArtifactFile(
            url = "https://caltech.box.com/shared/static/1uwg8rjg2wx7y0vp8j9kv2d44y3fajyk.csv",
            filename = "AMF_US-MOz_FLUXNET_FULLSET_HH_2005.csv",
        ),],
    )
    return AW.get_data_folder(ozark_fluxnet_dataset)
end

# LAI spatial data from MODIS; used for Ozark FLUXNET experiment
# TODO add citation
function ozark_lai_dataset_path()
    ozark_lai_dataset = AW.ArtifactWrapper(
        @__DIR__,
        "modis_data",
        AW.ArtifactFile[AW.ArtifactFile(
            url = "https://caltech.box.com/shared/static/yfqj0yqkx8yps7ydltsaixuy99083pon.csv",
            filename = "Ozark_MODIS_LAI_2005.csv",
        ),],
    )
    return AW.get_data_folder(ozark_lai_dataset)
end

# Leaf water potential (LWP) data; used for Ozark experiment
# TODO add citation
function ozark_lwp_dataset_path()
    ozark_lwp_dataset = AW.ArtifactWrapper(
        @__DIR__,
        "lwp_pallardy_etal2018",
        AW.ArtifactFile[AW.ArtifactFile(
            url = "https://caltech.box.com/shared/static/d2nbhezw1q99vslnh5qfwnqrnp3p4edo.csv",
            filename = "MOFLUX_PredawnLeafWaterPotential_2020_20210125.csv",
        ),],
    )
    return AW.get_data_folder(ozark_lwp_dataset)
end

# Soil water retention curve data; used for Ozark FLUXNET experiment
# [Full citation:] Wood, Jeffrey D, Gu, Lianhong, Hanson, Paul J, Frankenberg, Christian, & Sack,Lawren. (2022).
#  Supporting biophysical data for "The ecosystem wilting point defines drought response and recovery of a Quercus-Carya forest"
#  [Data set]. Zenodo. https://doi.org/10.5281/zenodo.7477879
function ozark_wrc_dataset_path()
    ozark_wrc_dataset = AW.ArtifactWrapper(
        @__DIR__,
        "ozark_swr_data",
        AW.ArtifactFile[AW.ArtifactFile(
            url = "https://caltech.box.com/shared/static/8aw3c3zygc7knw94py79826a9ai8hvkb.csv",
            filename = "MOFLUX_SWRC_data.csv",
        ),],
    )
    return AW.get_data_folder(ozark_wrc_dataset)
end

# Data from AmeriFlux Vaira FLUXNET site; used for Vaira FLUXNET experiment
# TODO add citation
function vaira_fluxnet_dataset_path()
    vaira_fluxnet_dataset = AW.ArtifactWrapper(
        @__DIR__,
        "ameriflux_data-US-Var",
        AW.ArtifactFile[AW.ArtifactFile(
            url = "https://caltech.box.com/shared/static/dx0p5idbsbrdebsda10t9pfv2lbdaz95.csv",
            filename = "AMF_US-Var_FLUXNET_FULLSET_HH_2003-2006.csv",
        ),],
    )
    return AW.get_data_folder(vaira_fluxnet_dataset)
end

# Evaporation data from Lehmann 2008; used for soil evaporation experiment
# [Full citation:] Lehmann et al, 2008. Fig. 8. Characteristic lengths affecting
#  evaporative drying of porous media. Phys. Rev. E 77, 056309.
function evap_dataset_path()
    evap_dataset = AW.ArtifactWrapper(
        @__DIR__,
        "lehmann2008_fig8_evaporation",
        AW.ArtifactFile[AW.ArtifactFile(
            url = "https://caltech.box.com/shared/static/cgppw3tx6zdz7h02yt28ri44g1j088ju.csv",
            filename = "lehmann2008_fig8_evaporation.csv",
        ),],
    )
    return AW.get_data_folder(evap_dataset)
end

# Data for clay soil type from Bonan 2019; used for soil Richards comparison experiment
# TODO citation
function bonan_clay_dataset_path()
    bonan_clay_dataset = AW.ArtifactWrapper(
        @__DIR__,
        "richards_clay",
        ArtifactFile[ArtifactFile(
            url = "https://caltech.box.com/shared/static/nk89znth59gcsdb65lnywnzjnuno3h6k.txt",
            filename = "clay_bonan_sp801_22323.txt",
        ),],
    )
    return AW.get_data_folder(bonan_clay_dataset)
end
# Data for sand soil type from Bonan 2019; used for soil Richards comparison experiment
# TODO citation
function bonan_sand_dataset_path()
    bonan_sand_dataset = ArtifactWrapper(
        @__DIR__,
        "richards_sand",
        ArtifactFile[ArtifactFile(
            url = "https://caltech.box.com/shared/static/2vk7bvyjah8xd5b7wxcqy72yfd2myjss.csv",
            filename = "sand_bonan_sp801.csv",
        ),],
    )
    return AW.get_data_folder(bonan_sand_dataset)
end

# Data for flux boundary condition case; used for Richards flux BC experiment
# TODO citation
function rre_flux_dataset_path()
    rre_flux_dataset = ArtifactWrapper(
        @__DIR__,
        "richards_flux_bc",
        ArtifactFile[ArtifactFile(
            url = "https://caltech.box.com/shared/static/bsfokpg0wvxoq04e8na0t3o0u6x5yw9n.csv",
            filename = "ref_soln_flux.csv",
        ),],
    )
    return AW.get_data_folder(rre_flux_dataset)
end

# Data for Dirichlet boundary condition case; used for Richards Dirichlet BC experiment
# TODO citation
function rre_dirichlet_dataset_path()
    rre_dirichlet_dataset = ArtifactWrapper(
        @__DIR__,
        "richards_dirichlet_bc",
        ArtifactFile[ArtifactFile(
            url = "https://caltech.box.com/shared/static/w6q30flbgj68lr0ncvoco10okrupmab1.csv",
            filename = "ref_soln_dirichlet.csv",
        ),],
    )
    return AW.get_data_folder(rre_dirichlet_dataset)
end


### TEST DATA
# Data for two stream radiation scheme; used for two stream unit tests
# TODO citation
function two_stream_dataset_path()
    two_stream_dataset = ArtifactWrapper(
        @__DIR__,
        "PySellersTwoStream Data",
        ArtifactFile[ArtifactFile(
            url = "https://caltech.box.com/shared/static/e7angzdnw18tmf8gctse5flsrkjsrlhx.csv",
            filename = "2_str_test_data.csv",
        ),],
    )
    return AW.get_data_folder(two_stream_dataset)
end


### BUCKET DATA
"""
    cesm2_albedo_dataset_path()

Triggers the download of the CESM2 land albedo dataset, if not
already downloaded, using Julia Artifacts, and returns the path to
this file.

This dataset contains monthly albedo data from 15/01/1850
to 15/12/2014.
"""
# TODO citation
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
# TODO citation
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

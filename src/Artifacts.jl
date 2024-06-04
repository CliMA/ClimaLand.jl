module Artifacts

import ClimaUtilities.ClimaArtifacts: @clima_artifact

import LazyArtifacts

using ArtifactWrappers

"""
    era5_land_forcing_data2021_path(; context)

Return the path to the folder that contains the ERA5 data.
"""
function era5_land_forcing_data2021_folder_path(; context = nothing)
    return @clima_artifact("era5_land_forcing_data2021", context)
end

"""
    soil_params_artifact_path(; context)

Return the path to the folder that contains the soil parameters.
"""
function soil_params_artifact_folder_path(; context = nothing)
    return @clima_artifact("soil_params_Gupta2020_2022", context)
end

"""
    experiment_fluxnet_data_path(
        site_ID,
        data_link;
        context = nothing,
    )

Return the path to the file that contains a year of fluxnet data
corresponding to a `site_ID` in the set (US-MOz, US-NR1, US-Ha1, US-Var).

Site publications and licenses:
US-Ha1: 

Urbanski, S., Barford, C., Wofsy, S., Kucharik, C., Pyle, E., Budney, J., McKain, K., Fitzjarrald, D., Czikowsky, M., Munger, J. W. 2007. Factors Controlling CO2 Exchange On Timescales From Hourly To Decadal At Harvard Forest, Journal Of Geophysical Research, 112:G2, . 

AmeriFlux FLUXNET: https://doi.org/10.17190/AMF/1871137
Citation: J. William Munger (2022), AmeriFlux FLUXNET-1F US-Ha1 Harvard Forest EMS Tower (HFR1), Ver. 3-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1871137

AmeriFlux CC-BY-4.0 Policy

US-MOz

Gu, L., Pallardy, S. G., Yang, B., Hosman, K. P., Mao, J., Ricciuto, D., Shi, X., Sun, Y. 2016. Testing a Land Model In Ecosystem Functional Space via a Comparison of Observed and Modeled Ecosystem Responses to Precipitation Regimes and Associated Stresses in a Central U.S. Forest, Journal Of Geophysical Research: Biogeosciences, 121:7, 1884-1902. 

AmeriFlux FLUXNET: https://doi.org/10.17190/AMF/1854370
Citation: Jeffrey Wood, Lianhong Gu (2021), AmeriFlux FLUXNET-1F US-MOz Missouri Ozark Site, Ver. 3-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1854370

AmeriFlux CC-BY-4.0 Policy

US-NR1

Burns, S. P., Blanken, P. D., Turnipseed, A. A., Hu, J., Monson, R. K. 2015. The Influence Of Warm-Season Precipitation On The Diel Cycle Of The Surface Energy Balance And Carbon Dioxide At A Colorado Subalpine Forest Site, Biogeosciences, 12:23, 7349-7377. 

AmeriFlux FLUXNET: https://doi.org/10.17190/AMF/1871141
Citation: Peter D. Blanken, Russel K. Monson, Sean P. Burns, David R. Bowling, Andrew A. Turnipseed (2022), AmeriFlux FLUXNET-1F US-NR1 Niwot Ridge Forest (LTER NWT1), Ver. 3-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1871141

AmeriFlux CC-BY-4.0 License

US-Var

AmeriFlux FLUXNET: https://doi.org/10.17190/AMF/1993904
Citation: Siyan Ma, Liukang Xu, Joseph Verfaillie, Dennis Baldocchi (2023), AmeriFlux FLUXNET-1F US-Var Vaira Ranch- Ione, Ver. 3-5, AmeriFlux AMP, (Dataset). https://doi.org/10.17190/AMF/1993904

AmeriFlux CC-BY-4.0 License
"""
function experiment_fluxnet_data_path(site_ID, data_link; context = nothing)
    @assert site_ID ∈ ("US-MOz", "US-Var", "US-NR1", "US-Ha1")
    dir = joinpath(@__DIR__, "../")
    af = ArtifactFile(
        url = data_link,
        filename = "AMF_$(site_ID)_FLUXNET_FULLSET.csv",
    )
    dataset =
        ArtifactWrapper(dir, "ameriflux_data_$(site_ID)", ArtifactFile[af])
    folder_path = get_data_folder(dataset)
    data_path = joinpath(folder_path, "AMF_$(site_ID)_FLUXNET_FULLSET.csv")
    return data_path
end

"""
    esm_snowmip_data_path(; site_ID = "cdp", context = nothing)

Returns the paths to the Col de Porte met data file and snow data file from the
ESM-SnowMIP data set.

Citation:
Menard, Cecile; Essery, Richard (2019): ESM-SnowMIP meteorological and evaluation datasets at ten reference sites (in situ and bias corrected reanalysis data) [dataset]. PANGAEA, https://doi.org/10.1594/PANGAEA.897575, Supplement to: Menard, Cecile; Essery, Richard; Barr, Alan; Bartlett, Paul; Derry, Jeff; Dumont, Marie; Fierz, Charles; Kim, Hyungjun; Kontu, Anna; Lejeune, Yves; Marks, Danny; Niwano, Masashi; Raleigh, Mark; Wang, Libo; Wever, Nander (2019): Meteorological and evaluation datasets for snow modelling at 10 reference sites: description of in situ and bias-corrected reanalysis data. Earth System Science Data, 11(2), 865-880, https://doi.org/10.5194/essd-11-865-2019

Creative Commons Attribution-NonCommercial 4.0 International (CC-BY-NC-4.0)

This dataset is soon to be in ClimaArtifacts and then we can update this function.
"""
function esm_snowmip_data_path(; site_ID = "cdp", context = nothing)
    dir = joinpath(@__DIR__, "../")
    met_af = ArtifactFile(
        url = "https://caltech.box.com/shared/static/q8ju0yy782j3f3jjrqowakpm5s3qf40v.nc",
        filename = "met_insitu_$(site_ID)_1994_2014.nc",
    )
    dataset = ArtifactWrapper(dir, "met drivers cdp", ArtifactFile[met_af])
    met_folder = get_data_folder(dataset)
    met_data_path = joinpath(met_folder, met_af.filename)

    snow_af = ArtifactFile(
        url = "https://caltech.box.com/shared/static/umwr6wqjhsz3zq6acht68m4d28a64ueb.nc",
        filename = "obs_insitu_$(site_ID)_1994_2014.nc",
    )
    snow_dataset = ArtifactWrapper(dir, "snow data cdp", ArtifactFile[snow_af])
    snow_folder = get_data_folder(snow_dataset)
    snow_data_path = joinpath(snow_folder, snow_af.filename)
    return met_data_path, snow_data_path
end

"""
    water_conservation_test_data_path(; context = nothing)

Returns the filepaths for data from two simulations of ClimaLand.Soil.RichardsModel;
these were carried out with a very small timestep with an explicit timestepper
and are used as ground truth for solutions using an implicit timestepper.

Experiment details are in `experiments/standalone/Soil/water_conservation.jl`.
"""
function water_conservation_test_data_path(; context = nothing)
    dir = joinpath(@__DIR__, "../")
    flux_dataset = ArtifactWrapper(
        dir,
        "richards_flux_bc_ref_soln",
        ArtifactFile[ArtifactFile(
            url = "https://caltech.box.com/shared/static/bsfokpg0wvxoq04e8na0t3o0u6x5yw9n.csv",
            filename = "ref_soln_flux.csv",
        ),],
    )
    flux_datapath = get_data_folder(flux_dataset)

    dirichlet_dataset = ArtifactWrapper(
        dir,
        "richards_dirichlet_bc_ref_soln",
        ArtifactFile[ArtifactFile(
            url = "https://caltech.box.com/shared/static/w6q30flbgj68lr0ncvoco10okrupmab1.csv",
            filename = "ref_soln_dirichlet.csv",
        ),],
    )
    dirichlet_datapath = get_data_folder(dirichlet_dataset)
    return flux_datapath, dirichlet_datapath
end

"""
    richards_eqn_bonan_data_path(; context = nothing)

Returns the file path for data created solving Richards equation
with G. Bonan's matlab code, found here:
https://github.com/gbonan/bonanmodeling/tree/master/sp_08_01

The data files correspond to the clay and sand data set described in
that code and in G. Bonan's book,
Climate Change and Terrestrial Ecosystem Modeling
DOI: https://doi.org/10.1017/9781107339217
Publisher: Cambridge University Press
Print publication year: 2019
"""
function richards_eqn_bonan_data_path(; context = nothing)
    dir = joinpath(@__DIR__, "../")
    bonan_clay_dataset = ArtifactWrapper(
        dir,
        "richards_clay_bonan_ref_soln",
        ArtifactFile[ArtifactFile(
            url = "https://caltech.box.com/shared/static/nk89znth59gcsdb65lnywnzjnuno3h6k.txt",
            filename = "clay_bonan_sp801_22323.txt",
        ),],
    )
    clay_datapath = joinpath(
        get_data_folder(bonan_clay_dataset),
        "clay_bonan_sp801_22323.txt",
    )

    bonan_sand_dataset = ArtifactWrapper(
        dir,
        "richards_sand_bonan_ref_soln",
        ArtifactFile[ArtifactFile(
            url = "https://caltech.box.com/shared/static/2vk7bvyjah8xd5b7wxcqy72yfd2myjss.csv",
            filename = "sand_bonan_sp801.csv",
        ),],
    )
    sand_datapath =
        joinpath(get_data_folder(bonan_sand_dataset), "sand_bonan_sp801.csv")
    return clay_datapath, sand_datapath
end

"""
    topmodel_data_path(; context = nothing)

Returns the path to the file which contains the necessary information for the TOPMODEL
runoff parameterization at 2.5 degrees resolution.

This file was created with 
https://github.com/CliMA/ClimaLand.jl/blob/main/src/standalone/Soil/Runoff/preprocess_topographic_index_simple.jl

using the data provided by
Marthews, T.R., Dadson, S.J., Lehner, B., Abele, S., Gedney, N. (2015). High-resolution global topographic index values. NERC Environmental Information Data Centre. (Dataset). https://doi.org/10.5285/6b0c4358-2bf3-4924-aa8f-793d468b92be

This resource is available under the Open Government Licence (OGL), and contains data supplied by Natural Environment Research Council.

This product, High-resolution global topographic index values, has been created with use of data from the HydroSHEDS database which is © World Wildlife Fund, Inc. (2006-2013) and has been used herein under license. The HydroSHEDS database and more information are available at http://www.hydrosheds.org.

Eventually, the script processing this data, and this data, will be added to ClimaArtifacts.
"""
function topmodel_data_path(; context = nothing)
    dir = joinpath(@__DIR__, "../")
    topmodel_dataset = ArtifactWrapper(
        dir,
        "processed_topographic_index 2.5 deg",
        ArtifactFile[ArtifactFile(
            url = "https://caltech.box.com/shared/static/dwa7g0uzhxd50a2z3thbx3a12n0r887s.nc",
            filename = "means_2.5_new.nc",
        ),],
    )
    path = joinpath(get_data_folder(topmodel_dataset), "means_2.5_new.nc")
    return path
end

end

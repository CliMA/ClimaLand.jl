module Artifacts

import ClimaUtilities.ClimaArtifacts: @clima_artifact

import LazyArtifacts


"""
    era5_land_forcing_data2021_path(; context)

Return the path to the folder that contains the ERA5 data
"""
function era5_land_forcing_data2021_folder_path(; context = nothing)
    return @clima_artifact("era5_land_forcing_data2021", context)
end

"""
    soil_params_artifact_path(; context)

Return the path to the folder that contains the soil parameters
"""
function soil_params_artifact_folder_path(; context = nothing)
    return @clima_artifact("soil_params_Gupta2020_2022", context)
end

end

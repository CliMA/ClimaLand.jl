include(joinpath(@__DIR__, "artifact_funcs.jl"))

# Trigger download if data doesn't exist locally
function trigger_download(lazy_download = true)
    # TUTORIAL DATA
    @info "Mizoguchi freezing front dataset path:`$(mizoguchi_dataset_path())`"

    # EXPERIMENT DATA
    @info "Ozark FLUXNET dataset path:`$(ozark_fluxnet_dataset_path())`"
    @info "Ozark LAI dataset path:`$(ozark_lai_dataset_path())`"
    @info "Ozark LWP dataset path:`$(ozark_lwp_dataset_path())`"
    @info "Ozark WRC dataset path:`$(ozark_wrc_dataset_path())`"
    @info "Vaira FLUXNET dataset path:`$(vaira_fluxnet_dataset_path())`"
    @info "Evaporation experiment dataset path:`$(evap_dataset_path())`"
    @info "Bonan clay dataset path:`$(bonan_clay_dataset_path())`"
    @info "Bonan sand dataset path:`$(bonan_sand_dataset_path())`"
    @info "Richards flux BC dataset path:`$(rre_flux_dataset_path())`"
    @info "Richards Dirichlet BC dataset path:`$(rre_dirichlet_dataset_path())`"

    # TEST DATA
    @info "Two stream dataset path:`$(two_stream_dataset_path())`"
    return nothing
end
trigger_download()

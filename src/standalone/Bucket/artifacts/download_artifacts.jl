include(joinpath(@__DIR__, "artifact_funcs.jl"))

# Trigger download if data doesn't exist locally
function trigger_download(lazy_download = true)
    @info "CESM2 temporal albedo dataset path:`$(cesm2_albedo_dataset_path())`"
    @info "Static albedo dataset path:`$(bareground_albedo_dataset_path())`"
    return nothing
end
trigger_download()

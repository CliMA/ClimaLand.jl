function get_path(site_ID)
    directory = "/groups/esm/ClimaArtifacts/artifacts/fluxnet_sites/"

    ret_str = directory * site_ID * ".csv"

    return ret_str
end

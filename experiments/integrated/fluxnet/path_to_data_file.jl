function get_path(site_ID)
    directories = readdir("/groups/esm/ClimaArtifacts/artifacts/fluxnet2015/")

    for directory in directories
        if occursin(site_ID, directory)
            parts = split(directory, r"FULLSET_", limit=2)

            parts[1] *= "FULLSET"
            ret_str = "/groups/esm/ClimaArtifacts/artifacts/fluxnet2015/" * directory * "/" * parts[1] * "_MM_" * parts[2] * ".csv"
            println(ret_str)
            return ret_str
        end
    end

end

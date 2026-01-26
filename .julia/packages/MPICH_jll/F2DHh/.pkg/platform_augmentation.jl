    using Base.BinaryPlatforms
    # Can't use Preferences since we might be running this very early with a non-existing Manifest
MPIPreferences_UUID = Base.UUID("3da0fdf6-3ccc-4f1b-acd9-58baa6c99267")
const preferences = Base.get_preferences(MPIPreferences_UUID)

# Keep logic in sync with MPIPreferences.jl
function augment_mpi!(platform)
    # Doesn't need to be `const` since we depend on MPIPreferences so we
    # invalidate the cache when it changes.
    # Note: MPIPreferences uses `Sys.iswindows()` without the `platform` argument.
    binary = get(preferences, "binary", Sys.iswindows(platform) ? "MicrosoftMPI_jll" : "MPICH_jll")

    abi = if binary == "system"
        let abi = get(preferences, "abi", nothing)
            if abi === nothing
                error("MPIPreferences: Inconsistent state detected, binary set to system, but no ABI set.")
            else
                abi
            end
        end
    elseif binary == "MicrosoftMPI_jll"
        "MicrosoftMPI"
    elseif binary == "MPICH_jll"
        "MPICH"
    elseif binary == "MPICH_CUDA_jll"
        "MPICH"
    elseif binary == "OpenMPI_jll"
        "OpenMPI"
    elseif binary == "MPItrampoline_jll"
        "MPItrampoline"
    else
        error("Unknown binary: $binary")
    end

    if !haskey(platform, "mpi")
        platform["mpi"] = abi
    end
    return platform
end

    augment_platform!(platform::Platform) = augment_mpi!(platform)

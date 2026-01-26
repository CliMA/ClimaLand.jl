module ScratchUsage
using Pkg, Scratch

# This function will create a bevy of spaces here
function touch_scratch()
    my_uuid = Base.PkgId(@__MODULE__).uuid
    my_version = Base.pkgversion(@__MODULE__)

    # Create an explicitly version-specific space
    private_space = get_scratch!(
        my_uuid,
        string(my_version.major, ".", my_version.minor, ".", my_version.patch),
    )
    touch(joinpath(private_space, string("ScratchUsage-", my_version)))

    # Create a space shared between all instances of the same major version,
    # using the `@get_scratch!` macro which automatically looks up the UUID
    major_space = @get_scratch!(string(my_version.major))
    touch(joinpath(major_space, string("ScratchUsage-", my_version)))

    # Create a global space that is not locked to this package at all
    # since this code is called from Pkg.test we need to pretend we are
    # at top level and not in the temporary Pkg.test project
    project = Base.ACTIVE_PROJECT[]
    Base.ACTIVE_PROJECT[] = "@v#.#"
    try
        global_space = get_scratch!("GlobalSpace")
        touch(joinpath(global_space, string("ScratchUsage-", my_version)))
    finally
        Base.ACTIVE_PROJECT[] = project
    end
end

end # module ScratchUsage

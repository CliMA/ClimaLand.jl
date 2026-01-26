# Things that are useful to know across platforms
if Sys.iswindows()
    const LIBPATH_env = "PATH"
    const LIBPATH_default = ""
    const pathsep = ';'
elseif Sys.isapple()
    const LIBPATH_env = "DYLD_FALLBACK_LIBRARY_PATH"
    const LIBPATH_default = "~/lib:/usr/local/lib:/lib:/usr/lib"
    const pathsep = ':'
else
    const LIBPATH_env = "LD_LIBRARY_PATH"
    const LIBPATH_default = ""
    const pathsep = ':'
end

function adjust_ENV!(env::Dict, PATH::String, LIBPATH::String, adjust_PATH::Bool, adjust_LIBPATH::Bool)
    if adjust_LIBPATH
        LIBPATH_base = get(env, LIBPATH_env, expanduser(LIBPATH_default))
        if !isempty(LIBPATH_base)
            env[LIBPATH_env] = string(LIBPATH, pathsep, LIBPATH_base)
        else
            env[LIBPATH_env] = LIBPATH
        end
    end
    if adjust_PATH && (LIBPATH_env != "PATH" || !adjust_LIBPATH)
        if adjust_PATH
            if !isempty(get(env, "PATH", ""))
                env["PATH"] = string(PATH, pathsep, env["PATH"])
            else
                env["PATH"] = PATH
            end
        end
    end
    return env
end

function withenv_executable_wrapper(f::Function,
                                    executable_path::String,
                                    PATH::String,
                                    LIBPATH::String,
                                    adjust_PATH::Bool,
                                    adjust_LIBPATH::Bool)
    env = Dict{String,String}(
        "PATH" => get(ENV, "PATH", ""),
        LIBPATH_env => get(ENV, LIBPATH_env, ""),
    )
    env = adjust_ENV!(env, PATH, LIBPATH, adjust_PATH, adjust_LIBPATH)
    withenv(env...) do
        f(executable_path)
    end
end

function dev_jll(src_name)
    # Grab `Pkg` module
    Pkg = first(filter(p-> p[1].name == "Pkg", Base.loaded_modules))[2]

    # First, `dev` out the package, but don't affect the current project
    mktempdir() do temp_env
        Pkg.activate(temp_env) do
            Pkg.develop("$(src_name)_jll")
            Pkg.instantiate()
        end
    end

    # Create the override directory by populating it with the artifact contents
    override_dir = joinpath(Pkg.devdir(), "$(src_name)_jll", "override")
    if !isdir(override_dir)
        artifacts_toml = joinpath(Pkg.devdir(), "$(src_name)_jll", "Artifacts.toml")
        art_hash = Pkg.Artifacts.artifact_hash(src_name, artifacts_toml)
        art_location = Pkg.Artifacts.artifact_path(art_hash)
        cp(art_location, override_dir)
    end
    # Force recompilation of that package, just in case it wasn't dev'ed before
    touch(joinpath(Pkg.devdir(), "$(src_name)_jll", "src", "$(src_name)_jll.jl"))
    @info("$(src_name)_jll dev'ed out to $(joinpath(Pkg.devdir(), "$(src_name)_jll")) with pre-populated override directory")
end


const JULIA_LIBDIRS = String[]
"""
    get_julia_libpaths()

Return the library paths that e.g. libjulia and such are stored in.
"""
function get_julia_libpaths()
    if isempty(JULIA_LIBDIRS)
        append!(JULIA_LIBDIRS, [joinpath(Sys.BINDIR::String, Base.LIBDIR, "julia"), joinpath(Sys.BINDIR::String, Base.LIBDIR)])
        # Windows needs to see the BINDIR as well
        @static if Sys.iswindows()
            push!(JULIA_LIBDIRS, Sys.BINDIR)
        end
    end
    return JULIA_LIBDIRS
end

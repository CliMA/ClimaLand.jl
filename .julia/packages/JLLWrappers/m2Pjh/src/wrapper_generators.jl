include("products/executable_generators.jl")
include("products/file_generators.jl")
include("products/library_generators.jl")

macro generate_wrapper_header(src_name)
    pkg_dir = dirname(dirname(String(__source__.file)))
    return esc(quote
        function find_artifact_dir()
            # We determine at compile-time whether our JLL package has been dev'ed and overridden
            @static if isdir(joinpath(dirname($(pkg_dir)), "override"))
                return joinpath(dirname($(pkg_dir)), "override")
            elseif @isdefined(augment_platform!) && VERSION >= v"1.6"
                $(Expr(:macrocall, Symbol("@artifact_str"), __source__, src_name, __module__.host_platform))
            else
                # We explicitly use `macrocall` here so that we can manually pass the `__source__`
                # argument, to avoid `@artifact_str` trying to lookup `Artifacts.toml` here.
                return $(Expr(:macrocall, Symbol("@artifact_str"), __source__, src_name))
            end
        end
        eager_mode() = nothing
        if ccall(:jl_generating_output, Cint, ()) == 1
            # to precompile these into Pkgimage
            Base.precompile(find_artifact_dir, ())
            Base.precompile(eager_mode, ())
        end
    end)
end


macro generate_init_header(dependencies...)
    deps_path_add = Expr[]
    eager_mode = Expr[]
    if !isempty(dependencies)
        for dep in dependencies
            depmod = getfield(__module__, dep)
            if isdefined(depmod, :PATH_list)
                push!(deps_path_add, :(append!(PATH_list, $(dep).PATH_list)))
            end
            if isdefined(depmod, :LIBPATH_list)
                push!(deps_path_add, :(append!(LIBPATH_list, $(dep).LIBPATH_list)))
            end
            if isdefined(depmod, :eager_mode)
                push!(eager_mode, :($(dep).eager_mode()))
            end
        end
    end

    return excat(
        # This either calls `@artifact_str()`, or returns a constant string if we're overridden.
        :(global artifact_dir = find_artifact_dir()::String),

        # Add `eager_mode` invocations on all our dependencies
        eager_mode...,

        # Initialize PATH_list and LIBPATH_list
        deps_path_add...,
    )
end


macro generate_init_footer()
    return esc(quote
        # Filter out duplicate and empty entries in our PATH and LIBPATH entries
        unique!(PATH_list)
        unique!(LIBPATH_list)
        PATH[] = join(PATH_list, $(pathsep))
        LIBPATH[] = join(vcat(LIBPATH_list, Base.invokelatest(JLLWrappers.get_julia_libpaths))::Vector{String}, $(pathsep))
    end)
end


"""
    emit_preference_path_load(pref_name, default_value)

On Julia 1.6+, emits a `load_preference()` call for the given preference name,
returning `nothing` if it is not loaded.  On Julia v1.5-, always returns `nothing`.
"""
function emit_preference_path_load(pref_name)
    # Can't use `Preferences.jl` on older Julias, just always use the default value in that case
    @static if VERSION < v"1.6.0-DEV"
        return quote
            nothing
        end
    else
        return quote
            @load_preference($(pref_name), nothing)
        end
    end
end


"""
    generate_imports(src_name)

We import the appropriate `Artifacts` module based on Julia version; Julia 1.3+ contains
an `Artifacts` module in `Pkg`, while Julia 1.6+ has it as a builtin standard library.
"""
function generate_imports(src_name)
    # We lie a bit in the registry that JLL packages are usable on Julia 1.0-1.2.
    # This is to allow packages that might want to support Julia 1.0 to get the
    # benefits of a JLL package on 1.3 (requiring them to declare a dependence on
    # the JLL package in their Project.toml) but engage in heroic hacks to do
    # something other than actually use a JLL package on 1.0-1.2.  By allowing
    # this package to be installed (but not loaded) on 1.0-1.2, we enable users
    # to avoid splitting their package versions into pre-1.3 and post-1.3 branches
    # if they are willing to engage in the kinds of hoop-jumping they might need
    # to in order to install binaries in a JLL-compatible way on 1.0-1.2. One
    # example of this hoop-jumping being to express a dependency on this JLL
    # package, then import it within a `VERSION >= v"1.3"` conditional, and use
    # the deprecated `build.jl` mechanism to download the binaries through e.g.
    # `BinaryProvider.jl`.  This should work well for the simplest packages, and
    # require greater and greater heroics for more and more complex packages.
    @static if VERSION < v"1.3.0-rc4"
        return quote
            error("Unable to use $($(src_name))_jll on Julia versions older than 1.3!")
        end
    elseif VERSION < v"1.6.0-DEV"
        # Use slow Pkg-based Artifacts
        return quote
            using Libdl, Pkg, Pkg.BinaryPlatforms, Pkg.Artifacts
            using Pkg.Artifacts: load_artifacts_toml, unpack_platform
            using Pkg.BinaryPlatforms: triplet, select_platform
            HostPlatform() = platform_key_abi()
        end
    else
        # Use fast stdlib-based Artifacts + Preferences
        return quote
            using Libdl, Artifacts, JLLWrappers.Preferences, Base.BinaryPlatforms
            using Artifacts: load_artifacts_toml, unpack_platform
            using Base.BinaryPlatforms: triplet, select_platform
        end
    end
end

"""
    generate_compiler_options(src_name)

Because JLL packages do not contain code that benefits much from compiler optimizations,
we disable them for a sizable boost in first load time.
"""
function generate_compiler_options(src_name)
    # Newer Julias have `@compiler_options` that can enable interpreted mode
    if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@compiler_options"))
        return quote
            Core.eval($(Symbol("$(src_name)_jll")), :(Base.Experimental.@compiler_options compile=min optimize=0 infer=false))
        end
    end

    # Older Julias only have `@optlevel`
    if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@optlevel"))
        return quote
            @eval Base.Experimental.@optlevel 0
        end
    end

    # If none of these are available, TOUGH BEANS.
    return nothing
end

"""
    generate_toplevel_definitions(src_name)

This method generates the toplevel definitions common to all JLL packages, such as
`is_available()`, the `PATH` and `LIBPATH` symbols, etc....
"""
function generate_toplevel_definitions(src_name, __source__)
    return quote
        """
            is_available()

        Return whether the artifact is available for the current platform.
        """
        function is_available end

        const PATH = Ref{String}("")
        const LIBPATH = Ref{String}("")
        # We put these inter-JLL-package API values here so that they are always defined, even if there
        # is no underlying wrapper held within this JLL package.
        const PATH_list = String[]
        const LIBPATH_list = String[]

        # We sub off to JLLWrappers' dev_jll, but avoid backedges
        dev_jll() = Base.invokelatest(JLLWrappers.dev_jll, $(src_name))
    end
end

"""
    generate_wrapper_load(src_name, pkg_uuid)

Because each platform could have completely disjoint products, we embed the information
for each in its own `.jl` file, named by triplet, and stored within the `src/wrappers`
directory.  This method generates the code to load `Artifacts.toml`, parse it for
artifacts, find the one that matches the host platform, then load the matching wrapper.
"""
function generate_wrapper_load(src_name, pkg_uuid, __source__)
    jll_name = "$(src_name)_jll"
    pkg_dir = dirname(String(__source__.file))

    function platform_parse_compat()
        @static if VERSION < v"1.6.0-DEV"
            return :(parse_wrapper_platform(x) = platform_key_abi(x))
        else
            # Use `tryparse` because the `Artifacts.toml` file cound include artifacts for
            # platforms/architectures not yet supported by the current version of Julia.
            return :(parse_wrapper_platform(x) = tryparse(Platform, x))
        end
    end

    return quote
        @static if $global_typeassert_available
            global best_wrapper::Union{Nothing,String}
        else
            global best_wrapper
        end
        # Load Artifacts.toml file and select best platform at compile-time, since this is
        # running at toplevel, and therefore will be run completely at compile-time.  We use
        # a `let` block here to avoid storing unnecessary data in our `.ji` files

        if @isdefined(augment_platform!)
            const host_platform = augment_platform!(HostPlatform())
        else
            const host_platform = nothing
        end
        best_wrapper = let
            artifacts_toml = joinpath($(pkg_dir), "..", "Artifacts.toml")
            valid_wrappers = Dict{Platform,String}()
            artifacts = load_artifacts_toml(artifacts_toml; pkg_uuid=$(pkg_uuid))[$(src_name)]

            # Helper function to parse triplets for us
            $(platform_parse_compat())
            function make_wrapper_dict(dir, x)
                # (1) make the Dict type inferrable
                # (2) avoid creation of Generators so that we don't have to compile
                #     package-specific Generator-closures
                d = Dict{Platform,String}()
                for f in x
                    platform = parse_wrapper_platform(basename(f)[1:end-3])
                    # `platform` could be `nothing` for example if the architecture isn't
                    # known to currently running Julia version.
                    if !isnothing(platform)
                        d[platform] = joinpath(dir, "wrappers", f)
                    end
                end
                return d
            end

            # If it's a Dict, that means this is an AnyPlatform artifact, act accordingly:
            if isa(artifacts, Dict)
                joinpath($(pkg_dir), "wrappers", "any.jl")
            else
                # Otherwise, it's a Vector, and we must select the best platform
                # First, find all wrappers on-disk, parse their platforms, and match:
                wrapper_files = String[]
                for x in readdir(joinpath($(pkg_dir), "wrappers"))
                    endswith(x, ".jl") && push!(wrapper_files, x)   # avoid creation of Generators...
                end
                wrappers = make_wrapper_dict($(pkg_dir), wrapper_files)
                for e in artifacts
                    platform = unpack_platform(e, $(jll_name), artifacts_toml)

                    # Filter platforms based on what wrappers we've generated on-disk.
                    # Because the wrapper file naming strategy relies upon BB's assumptions of triplet
                    # parsing, it can be somewhat fragile.  So we have loaded all the wrapper filenames
                    # in, parsed them, and are now using `select_platform()` to match them, which is
                    # much more robust.  This step avoids a disconnect between what is recorded in the
                    # `Artifacts.toml` file and what wrappers are available on-disk.
                    wrapper_file = select_platform(wrappers, platform)
                    if wrapper_file !== nothing
                        valid_wrappers[platform] = wrapper_file
                    end
                end

                # From the available options, choose the best wrapper script
                # The two argument `select_platform` is notably slower, so micro-optimize this by
                # only calling it when necessary.
                if host_platform !== nothing
                    select_platform(valid_wrappers, host_platform)
                else
                    select_platform(valid_wrappers)
                end
            end
        end

        # Load in the wrapper, if it's not `nothing`!
        if best_wrapper === nothing
            @debug(string("Unable to load ", $(src_name), "; unsupported platform ", host_platform === nothing ? triplet(HostPlatform()) : triplet(host_platform)))
            is_available() = false
        else
            Base.include($(Symbol("$(src_name)_jll")), best_wrapper)
            is_available() = true
        end
    end
end

macro generate_main_file_header(src_name)
    return excat(
        # Declare this module as interpreted
        generate_compiler_options(src_name),
        # import Artifacts module
        generate_imports(src_name),
        global_typeassert_available ? :(global artifact_dir::String) : :()
    )
end

"""
    @generate_main_file(src_name, pkg_uuid)

Generate the main JLL file to Import `Artifacts`, set compiler options, perform forward
declarations for JLL-internal APIs, load the wrappers, etc...
"""
macro generate_main_file(src_name, pkg_uuid)
    return excat(
        # `is_available()` forward declaration, `PATH_list` and `LIBPATH_list` forward definitions,
        # `find_artifact_dir()` definition.
        generate_toplevel_definitions(src_name, __source__),
        # Select and load the best wrapper file
        generate_wrapper_load(src_name, pkg_uuid, __source__),
    )
end

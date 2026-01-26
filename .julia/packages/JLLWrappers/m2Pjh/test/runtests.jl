using JLLWrappers
using Pkg
using Test
using Downloads
using Artifacts

@static if VERSION >= v"1.6.0-DEV"
    using Preferences
end

# Older versions of Julia need to test versions without `aarch64-apple-darwin` mappings,
# because older Julia was too picky about what kinds of platforms you're allowed to create.
if VERSION <= v"1.6.0-DEV"
    test_pkgs = [
        ("HelloWorldC", v"1.3.0+noapple"),
        ("libxls", v"1.6.2+noapple"),
        ("Vulkan_Headers", v"1.3.243+1"),
    ]
else
    test_pkgs = [
        ("HelloWorldC", v"1.4.0+0"),
        ("LAMMPS", v"2.4.0+0"),
        ("libxls", v"1.6.2+0"),
        ("Vulkan_Headers", v"1.3.243+1"),
    ]
end

mktempdir() do dir; Pkg.activate(dir) do
    for (name, version) in test_pkgs
        path = joinpath(@__DIR__, "$(name)_jll")
        if !isdir(joinpath(path, "src"))
            rm(path; recursive=true, force=true)

            @info("Downloading test package", name, version)
            filename = "$(name)-v$(version).tar.gz"
            url = "https://github.com/JuliaBinaryWrappers/$(name)_jll.jl/archive/$(filename)"
            mkpath(path)
            Downloads.download(url, joinpath(path, filename))
            run(Cmd(`tar --strip-components=1 -zxf $(filename)`; dir=path))
        end
    end
end; end

module TestJLL end
@testset "JLLWrappers.jl" begin
    mktempdir() do dir
        Pkg.activate(dir)

        # actually use the development version of JLLWrappers
        Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))

        # Manually calculate the artifacts path for `libxlsreader` so that we don't actually load the JLL,
        # because we're about to overwrite it
        libxlsreader_artifact_dir = Artifacts.artifact_path(Base.SHA1(Artifacts.artifact_meta("libxls", "./libxls_jll/Artifacts.toml")["git-tree-sha1"]))
        libxlsreader_unversioned_path = joinpath(
            libxlsreader_artifact_dir,
            Sys.iswindows() ? "bin" : "lib",
            string("libxlsreader.", Pkg.BinaryPlatforms.platform_dlext()),
        )

        # Prepare some overrides for various products
        @static if VERSION >= v"1.6.0-DEV"
            set_preferences!(joinpath(dir, "LocalPreferences.toml"), "Vulkan_Headers_jll", "vulkan_hpp_path" => "foo")
            set_preferences!(joinpath(dir, "LocalPreferences.toml"), "libxls_jll", "libxlsreader_path" => libxlsreader_unversioned_path)
            set_preferences!(joinpath(dir, "LocalPreferences.toml"), "HelloWorldC_jll", "hello_world_doppelganger_path" => "foo")
        end

        # Package with a FileProduct
        Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "Vulkan_Headers_jll")))
        @eval TestJLL using Vulkan_Headers_jll
        @test @eval TestJLL Vulkan_Headers_jll.is_available()
        @test isfile(@eval TestJLL vk_xml)
        @test isfile(@eval TestJLL Vulkan_Headers_jll.vk_xml_path)
        @test isfile(@eval TestJLL Vulkan_Headers_jll.get_vk_xml_path())
        @static if VERSION >= v"1.6.0-DEV"
            @test !isfile(@eval TestJLL vulkan_hpp)
            @test !isfile(@eval TestJLL Vulkan_Headers_jll.vulkan_hpp_path)
            @test !isfile(@eval TestJLL Vulkan_Headers_jll.get_vulkan_hpp_path())
        else
            # If we don't have Preferences support, we should still find these.
            @test isfile(@eval TestJLL vulkan_hpp)
            @test isfile(@eval TestJLL Vulkan_Headers_jll.vulkan_hpp_path)
            @test isfile(@eval TestJLL Vulkan_Headers_jll.get_vulkan_hpp_path())
        end
        @test isdir(@eval TestJLL Vulkan_Headers_jll.artifact_dir)
        @test isempty(@eval TestJLL Vulkan_Headers_jll.PATH[])
        @test occursin(Sys.BINDIR, @eval TestJLL Vulkan_Headers_jll.LIBPATH[])

        # Package with an ExecutableProduct
        Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "HelloWorldC_jll")))
        @eval TestJLL using HelloWorldC_jll
        @test @eval TestJLL HelloWorldC_jll.is_available()
        @test isfile(@eval TestJLL HelloWorldC_jll.hello_world_path)
        @test isfile(@eval TestJLL HelloWorldC_jll.get_hello_world_path())
        @test startswith(basename(@eval TestJLL HelloWorldC_jll.get_hello_world_path()), "hello_world")
        @test isdir(@eval TestJLL HelloWorldC_jll.artifact_dir)
        @test !isempty(@eval TestJLL HelloWorldC_jll.PATH_list)
        @test occursin(Sys.BINDIR, @eval TestJLL HelloWorldC_jll.LIBPATH[])

        @static if VERSION >= v"1.6.0-DEV"
            @test !isfile(@eval TestJLL HelloWorldC_jll.hello_world_doppelganger_path)
            @test !isfile(@eval TestJLL HelloWorldC_jll.get_hello_world_doppelganger_path())
            @test @eval TestJLL HelloWorldC_jll.get_hello_world_doppelganger_path() == "foo"
        else
            @test isfile(@eval TestJLL HelloWorldC_jll.hello_world_doppelganger_path)
            @test isfile(@eval TestJLL HelloWorldC_jll.get_hello_world_doppelganger_path())
            @test @eval TestJLL startswith(basename(HelloWorldC_jll.get_hello_world_doppelganger_path()), "hello_world")
        end

        # Package with a LibraryProduct
        Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "libxls_jll")))
        # Windows doesn't have both, so we have to make sure it exists:
        if Sys.iswindows() && !isfile(joinpath(libxlsreader_artifact_dir, "bin", "libxlsreader.dll"))
            cp(
                joinpath(libxlsreader_artifact_dir, "bin", "libxlsreader-8.dll"),
                joinpath(libxlsreader_artifact_dir, "bin", "libxlsreader.dll")
            )
            chmod(joinpath(libxlsreader_artifact_dir, "bin", "libxlsreader.dll"), 0o755)
        end
        @eval TestJLL using libxls_jll
        @test @eval TestJLL libxls_jll.is_available()
        @test isfile(@eval TestJLL libxls_jll.libxlsreader_path)
        @test isfile(@eval TestJLL libxls_jll.get_libxlsreader_path())
        @test isdir(@eval TestJLL libxls_jll.artifact_dir)
        @test isempty(@eval TestJLL libxls_jll.PATH[])
        @test occursin(Sys.BINDIR, @eval TestJLL libxls_jll.LIBPATH[])
        libxlsreader_path = @eval TestJLL libxls_jll.get_libxlsreader_path()
        @static if VERSION >= v"1.6.0-DEV"
            @test endswith(libxlsreader_path, libxlsreader_unversioned_path)
        end

        # Issue #20
        if Sys.iswindows()
            @test Sys.BINDIR ∈ JLLWrappers.get_julia_libpaths()
        end

        # Package with an augment_platform!
        @static if VERSION >= v"1.6.0-DEV"
            Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "LAMMPS_jll")))
            if VERSION < v"1.7"
                # On 1.6 we lazyily download the missing artifact, causing output to stderr
                @eval TestJLL using LAMMPS_jll
            else
                @eval TestJLL using LAMMPS_jll
            end
            if @eval TestJLL LAMMPS_jll.is_available()
                mpi_abi = @eval TestJLL LAMMPS_jll.augment_platform!(LAMMPS_jll.HostPlatform())["mpi"]
                @test mpi_abi == lowercase(Sys.iswindows() ? "MicrosoftMPI" : "mpich")
                @test isfile(@eval TestJLL LAMMPS_jll.liblammps_path)
                @test isfile(@eval TestJLL LAMMPS_jll.get_liblammps_path())
                @test isdir(@eval TestJLL LAMMPS_jll.artifact_dir)
                @test occursin(Sys.BINDIR, @eval TestJLL LAMMPS_jll.LIBPATH[])
            end
        end
    end
end

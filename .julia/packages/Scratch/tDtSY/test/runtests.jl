using Scratch, Test, Dates, Pkg, Pkg.TOML
include("utils.jl")

# Set to true for verbose Pkg output
const verbose = false
global const pkgio = verbose ? stderr : devnull
@testset "Scratch Space Basics" begin
    # Run everything in a separate depot, so that we can test GC'ing and whatnot
    temp_pkg_dir() do project_dir
        # Create a global scratch space, ensure it exists and is writable
        dir = get_scratch!("test")
        @test isdir(dir)
        @test startswith(dir, scratch_dir())
        touch(joinpath(dir, "foo"))
        @test readdir(dir) == ["foo"]

        # Test that this created a `scratch_usage.toml` file, and that accessing it
        # again does not increase the size of the scratch_usage.toml file, since we
        # only mark usage once every so often per julia session.
        usage_path = Scratch.usage_toml()
        @test isfile(usage_path)
        size = filesize(usage_path)
        dir = get_scratch!("test")
        @test size == filesize(usage_path)

        # But accessing it from a new Julia instance WILL increase its size:
        code = "import Scratch; Scratch.get_scratch!(\"test\")"
        run(setenv(
            `$(Base.julia_cmd()) --project=$(dirname(@__DIR__)) -e $code`,
            "JULIA_DEPOT_PATH" => first(Base.DEPOT_PATH),
        ))
        @test size < filesize(usage_path)

        # Delete the scratch space, ensure it's gone.
        delete_scratch!("test")
        @test !isdir(dir)

        # Key verification
        @test isdir(get_scratch!("abcABC123._-"))
        @test_throws ArgumentError get_scratch!("")
        @test_throws ArgumentError get_scratch!("hello/world")
        @test_throws ArgumentError get_scratch!("hello\\world")
    end
end

@testset "Scratch Space Namespacing" begin
    temp_pkg_dir() do project_dir
        # ScratchUsage UUID
        su_uuid = "93485645-17f1-6f3b-45bc-419db53815ea"
        # The UUID that gets used when no good UUIDs are available
        global_uuid = string(Base.UUID(UInt128(0)))

        # Touch the spaces of a ScratchUsage v1.0.0
        install_test_ScratchUsage(project_dir, v"1.0.0")

        # Ensure that the files were created for v1.0.0
        @test isfile(scratch_dir(su_uuid, "1.0.0", "ScratchUsage-1.0.0"))
        @test length(readdir(scratch_dir(su_uuid, "1.0.0"))) == 1
        @test isfile(scratch_dir(su_uuid, "1", "ScratchUsage-1.0.0"))
        @test length(readdir(scratch_dir(su_uuid, "1"))) == 1
        @test isfile(scratch_dir(global_uuid, "GlobalSpace", "ScratchUsage-1.0.0"))
        @test length(readdir(scratch_dir(global_uuid, "GlobalSpace"))) == 1

        # Next, do the same but for more versions
        install_test_ScratchUsage(project_dir, v"1.1.0")
        install_test_ScratchUsage(project_dir, v"2.0.0")

        # Check the spaces were shared when they should have been, and not when they shouldn't
        @test isfile(scratch_dir(su_uuid, "1.0.0", "ScratchUsage-1.0.0"))
        @test length(readdir(scratch_dir(su_uuid, "1.0.0"))) == 1
        @test isfile(scratch_dir(su_uuid, "1.1.0", "ScratchUsage-1.1.0"))
        @test length(readdir(scratch_dir(su_uuid, "1.1.0"))) == 1
        @test isfile(scratch_dir(su_uuid, "2.0.0", "ScratchUsage-2.0.0"))
        @test length(readdir(scratch_dir(su_uuid, "2.0.0"))) == 1
        @test isfile(scratch_dir(su_uuid, "1", "ScratchUsage-1.0.0"))
        @test isfile(scratch_dir(su_uuid, "1", "ScratchUsage-1.1.0"))
        @test length(readdir(scratch_dir(su_uuid, "1"))) == 2
        @test isfile(scratch_dir(su_uuid, "2", "ScratchUsage-2.0.0"))
        @test length(readdir(scratch_dir(su_uuid, "2"))) == 1
        @test isfile(scratch_dir(global_uuid, "GlobalSpace", "ScratchUsage-1.0.0"))
        @test isfile(scratch_dir(global_uuid, "GlobalSpace", "ScratchUsage-1.1.0"))
        @test isfile(scratch_dir(global_uuid, "GlobalSpace", "ScratchUsage-2.0.0"))
        @test length(readdir(scratch_dir(global_uuid, "GlobalSpace"))) == 3

        clear_scratchspaces!(Base.UUID(su_uuid))
        @test !isdir(scratch_dir(su_uuid))

        # UUID lookup for active project when running in Main
        ## Project.toml without UUID
        project = temp_project_file()
        with_active_project(project) do
            @test (@__MODULE__) == Main
            ## Project.toml without UUID
            path = @get_scratch!("project-no-uuid")
            @test isdir(path)
            @test path == scratch_dir(global_uuid, "project-no-uuid")
            @test path === get_scratch!(@__MODULE__, "project-no-uuid")
        end
        ## Project.toml with UUID
        project_uuid = Base.UUID("69386cca-e009-4a96-a0ae-829213699cfc")
        project = temp_project_file(project_uuid)
        with_active_project(project) do
            path = @get_scratch!("project-uuid")
            @test isdir(path)
            @test path == scratch_dir(string(project_uuid), "project-uuid")
            @test path === get_scratch!(@__MODULE__, "project-uuid")

            ## Issue #29; when the Project.toml is not accessible, we should not throw an error
            if Sys.isunix()
                chmod(dirname(project), 0o000)
                Scratch.track_scratch_access(project_uuid, "not a real path")
                chmod(dirname(project), 0o755)
            end
        end # do

        # Cross-package scratch usage: Test that the scratch space is namespaced
        # to the other package, but tracked to me when it comes to lifecycling
        # the other package uuid, but the scratchspace is
        other_uuid = Base.UUID("6dc9890c-246d-42f7-b07c-3ce39ca50d56")
        ## Activate test-project from above again
        with_active_project(project) do
            path = get_scratch!(other_uuid, "hello-there", project_uuid)
            ## Test that the path is namespaced to other_uuid, but lifecycled with me
            @test path == scratch_dir(string(other_uuid), "hello-there")
            usage_path = usage_path = Scratch.usage_toml()
            usage = Pkg.TOML.parsefile(usage_path)
            @test project ∈ usage[path][1]["parent_projects"]
        end

        # Internal scratch access tracking
        empty!(Scratch.scratch_access_timers)
        usage_path = usage_path = Scratch.usage_toml()
        project2 = temp_project_file(other_uuid)
        with_active_project(project) do
            ## This should track project_uuid as the user
            path = get_scratch!(project_uuid, "general-kenobi")
            @test project_uuid ∈ first.(keys(Scratch.scratch_access_timers))
            usage = Pkg.TOML.parsefile(usage_path)
            @test any(project ∈ record["parent_projects"] for record in usage[path])
        end
        with_active_project(project2) do
            ## Reaching for the same space as another owner should (i) track this
            ## usage and (ii) write to scratch_usage.toml since this belongs to
            ## another project file
            path = get_scratch!(project_uuid, "general-kenobi", other_uuid)
            @test other_uuid ∈ first.(keys(Scratch.scratch_access_timers))
            usage = Pkg.TOML.parsefile(usage_path)
            @test any(project2 ∈ record["parent_projects"] for record in usage[path])
        end
        ## Deleting the path should remove both UUIDs from the timers
        delete_scratch!(project_uuid, "general-kenobi")
        @test project_uuid ∉ first.(keys(Scratch.scratch_access_timers))
        @test other_uuid ∉ first.(keys(Scratch.scratch_access_timers))
    end
end

@testset "Scratch Space Lifecycling" begin
    temp_pkg_dir() do project_dir
        # First, install ScratchUsage
        su_uuid = "93485645-17f1-6f3b-45bc-419db53815ea"
        global_uuid = string(Base.UUID(UInt128(0)))
        install_test_ScratchUsage(project_dir, v"1.0.0")

        # Ensure that a few files were created
        @test isfile(scratch_dir(su_uuid, "1.0.0", "ScratchUsage-1.0.0"))
        @test length(readdir(scratch_dir(su_uuid, "1.0.0"))) == 1
        @test isfile(scratch_dir(global_uuid, "GlobalSpace", "ScratchUsage-1.0.0"))
        @test length(readdir(scratch_dir(global_uuid, "GlobalSpace"))) == 1

        # Test that a gc() doesn't remove anything, and that there is no orphanage
        Pkg.gc(; io=pkgio)
        orphaned_path = joinpath(first(Base.DEPOT_PATH), "logs", "orphaned.toml")
        @test isfile(scratch_dir(su_uuid, "1.0.0", "ScratchUsage-1.0.0"))
        @test isfile(scratch_dir(global_uuid, "GlobalSpace", "ScratchUsage-1.0.0"))
        @test !isfile(orphaned_path) || filesize(orphaned_path) == 0

        # Remove ScrachUsage, which causes the package (but not the scratch dirs)
        # to move to the orphanage
        Pkg.rm("ScratchUsage"; io=pkgio)
        rm(joinpath(project_dir, "ScratchUsage"); force=true, recursive=true)
        Pkg.gc(; io=pkgio)

        @test isfile(scratch_dir(su_uuid, "1.0.0", "ScratchUsage-1.0.0"))
        @test isfile(scratch_dir(global_uuid, "GlobalSpace", "ScratchUsage-1.0.0"))
        @test isfile(orphaned_path)
        orphanage = Pkg.TOML.parse(String(read(orphaned_path)))
        @test haskey(orphanage, scratch_dir(su_uuid, "1.0.0"))
        @test haskey(orphanage, scratch_dir(su_uuid, "1"))
        @test !haskey(orphanage, scratch_dir(global_uuid, "GlobalSpace"))

        # Run a GC, forcing collection to ensure that everything in the SpaceUsage
        # namespace gets removed (but still appears in the orphanage)
        sleep(0.2)
        Pkg.gc(;collect_delay=Millisecond(100), io=pkgio)
        @test !isdir(scratch_dir(su_uuid))
        @test isdir(scratch_dir(global_uuid, "GlobalSpace"))
        orphanage = Pkg.TOML.parse(String(read(orphaned_path)))
        @test haskey(orphanage, scratch_dir(su_uuid, "1.0.0"))
        @test haskey(orphanage, scratch_dir(su_uuid, "1"))
        @test !haskey(orphanage, scratch_dir(global_uuid, "GlobalSpace"))
    end
end

@testset "Scratch Space Advanced Keywords" begin
    mktempdir() do depot_path
        # Test that targeting via `depot_path` works
        path = get_scratch!("foo"; depot_path)
        @test startswith(path, depot_path)

        # Test that this wrote the scratch access in the right place
        usage = TOML.parsefile(Scratch.usage_toml(depot_path))
        @test length(usage[path]) == 1

        # Test that the first call gets swalled due to the time gate,
        # but the second makes it through because we've disabled it:
        get_scratch!("foo"; depot_path)
        usage = TOML.parsefile(Scratch.usage_toml(depot_path))
        @test length(usage[path]) == 1

        get_scratch!("foo"; depot_path, time_gate=Second(-1))
        usage = TOML.parsefile(Scratch.usage_toml(depot_path))
        @test length(usage[path]) == 2

        # Test that the third call gets through as well, because
        # we lie and say the time is beyond the time gate:
        get_scratch!("foo"; depot_path, time_gate=Second(10), curr_time = Dates.now() + Second(11))
        usage = TOML.parsefile(Scratch.usage_toml(depot_path))
        @test length(usage[path]) == 3
    end
end


# Run a test using JET to do some static analysis for us
if Base.VERSION >= v"1.10"
    using JET
    @testset "test_package" begin
        test_package("Scratch")
    end
end

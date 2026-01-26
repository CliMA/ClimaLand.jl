import ClimaUtilities.OutputPathGenerator:
    generate_output_path,
    RemovePreexistingStyle,
    ActiveLinkStyle,
    detect_restart_file
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import Base: rm
using Test

const context = ClimaComms.context()
ClimaComms.init(context)

let_filesystem_catch_up() = context isa ClimaComms.MPICommsContext && sleep(0.2)

@testset "RemovePrexistingStyle" begin
    # Test empty output_path
    @test_throws ErrorException generate_output_path(
        "",
        style = RemovePreexistingStyle(),
    )

    base_output_path = ClimaComms.iamroot(context) ? mktempdir(pwd()) : ""
    base_output_path = ClimaComms.bcast(context, base_output_path)
    ClimaComms.barrier(context)
    let_filesystem_catch_up()

    # Folder does not yet exist
    output_path = joinpath(base_output_path, "dormouse")
    @test output_path == generate_output_path(
        output_path,
        context = context,
        style = RemovePreexistingStyle(),
    )

    # Check that it exists now
    @test isdir(output_path)

    if ClimaComms.iamroot(context)
        # Now the folder exists, let us add a file there
        open(joinpath(output_path, "something"), "w") do file
            write(file, "Something")
        end
    end
    ClimaComms.barrier(context)
    let_filesystem_catch_up()
    @test isfile(joinpath(output_path, "something"))

    # Check that the file got removed
    @test output_path == generate_output_path(
        output_path,
        context = context,
        style = RemovePreexistingStyle(),
    )
    let_filesystem_catch_up()

    @test !isfile(joinpath(output_path, "something"))

    if ClimaComms.iamroot(context)
        Base.rm(base_output_path, force = true, recursive = true)
    end

    ClimaComms.barrier(context)

    let_filesystem_catch_up()
end

@testset "ActiveLinkStyle" begin
    # Test empty output_path
    @test_throws ErrorException generate_output_path("")

    base_output_path = ClimaComms.iamroot(context) ? mktempdir(pwd()) : ""
    base_output_path = ClimaComms.bcast(context, base_output_path)
    ClimaComms.barrier(context)
    let_filesystem_catch_up()

    # Folder does not yet exist
    output_path = joinpath(base_output_path, "dormouse")

    output_link = joinpath(output_path, "output_active")

    expected_output = joinpath(output_path, "output_0000")

    @test expected_output ==
          generate_output_path(output_path, context = context)

    # Check that it exists now
    @test isdir(output_path)

    # Check link output_active was created
    @test islink(output_link)

    # # Check link points to folder
    @test readlink(output_link) == "output_0000"

    # Now the folder exists, let us see if the rotation works
    expected_output = joinpath(output_path, "output_0001")

    @test expected_output ==
          generate_output_path(output_path, context = context)

    # Check folder was created
    @test isdir(joinpath(output_path, "output_0001"))

    # Check link points to new folder
    @test readlink(output_link) == "output_0001"

    # Now let us check something wrong

    # Missing link and existing output_ folders
    if ClimaComms.iamroot(context)
        rm(expected_output)
    end
    ClimaComms.barrier(context)
    let_filesystem_catch_up()

    generate_output_path(output_path, context = context)
    @test readlink(output_link) == "output_0002"

    ClimaComms.barrier(context)
    let_filesystem_catch_up()

    # Remove link, we are going to create a new one
    if ClimaComms.iamroot(context)
        Base.rm(output_link)
    end

    ClimaComms.barrier(context)
    let_filesystem_catch_up()

    # Wrong link
    if ClimaComms.iamroot(context)
        wrong_dir = joinpath(output_path, "wrong")
        mkdir(wrong_dir)
        symlink(wrong_dir, output_link, dir_target = true)
    end
    ClimaComms.barrier(context)
    @test_throws ErrorException generate_output_path(
        output_path,
        context = context,
    )

    if ClimaComms.iamroot(context)
        Base.rm(base_output_path, force = true, recursive = true)
    end
    ClimaComms.barrier(context)
    let_filesystem_catch_up()
end

@testset "detect_restart_file" begin
    # Test with non-ActiveLinkStyle
    @test_throws ErrorException detect_restart_file(
        ".",
        style = RemovePreexistingStyle(),
    )

    # Test with non-existent base directory
    base_output_dir = "non_existent_dir"
    restart_file = detect_restart_file(base_output_dir)
    @test isnothing(restart_file)

    # Test with empty base directory
    mktempdir(pwd()) do base_output_dir
        restart_file = detect_restart_file(base_output_dir)
        @test isnothing(restart_file)
    end

    # Test with a single output directory and no restart files
    mktempdir(pwd()) do base_output_dir
        mkdir(joinpath(base_output_dir, "output_0001"))
        restart_file = detect_restart_file(base_output_dir)
        @test isnothing(restart_file)
    end

    # Test two output directories and no restart files
    mktempdir(pwd()) do base_output_dir
        mkdir(joinpath(base_output_dir, "output_0001"))
        mkdir(joinpath(base_output_dir, "output_0002"))
        restart_file = detect_restart_file(base_output_dir)
        @test isnothing(restart_file)
    end

    # Test with a single output directory and a single restart file
    mktempdir(pwd()) do base_output_dir
        output_dir = joinpath(base_output_dir, "output_0001")
        mkdir(output_dir)
        touch(joinpath(output_dir, "day0001.1234.hdf5"))
        restart_file = detect_restart_file(base_output_dir)
        @test restart_file == joinpath(output_dir, "day0001.1234.hdf5")
    end

    # Test with multiple output directories and restart files
    mktempdir(pwd()) do base_output_dir
        output_dir1 = joinpath(base_output_dir, "output_0001")
        output_dir2 = joinpath(base_output_dir, "output_0002")
        mkdir(output_dir1)
        mkdir(output_dir2)
        touch(joinpath(output_dir1, "day0001.1234.hdf5"))
        sleep(0.1)
        touch(joinpath(output_dir2, "day0002.5678.hdf5"))
        restart_file = detect_restart_file(base_output_dir)
        @test restart_file == joinpath(output_dir2, "day0002.5678.hdf5")
    end

    # Test with multiple output directories and restart files only in the first
    mktempdir(pwd()) do base_output_dir
        output_dir1 = joinpath(base_output_dir, "output_0001")
        output_dir2 = joinpath(base_output_dir, "output_0002")
        mkdir(output_dir1)
        mkdir(output_dir2)
        touch(joinpath(output_dir1, "day0001.1234.hdf5"))
        restart_file = detect_restart_file(base_output_dir)
        @test restart_file == joinpath(output_dir1, "day0001.1234.hdf5")
    end

    # Test with a custom restart file regular expression
    mktempdir(pwd()) do base_output_dir
        output_dir = joinpath(base_output_dir, "output_0001")
        mkdir(output_dir)
        touch(joinpath(output_dir, "restart_0001.h5"))
        restart_file = detect_restart_file(
            ActiveLinkStyle(),
            base_output_dir;
            restart_file_rx = r"restart_\d+\.h5",
        )
        @test restart_file == joinpath(output_dir, "restart_0001.h5")
    end

    # Test with a custom sorting function (e.g., sort by reverse modification time)
    mktempdir(pwd()) do base_output_dir
        output_dir = joinpath(base_output_dir, "output_0001")
        mkdir(output_dir)
        touch(joinpath(output_dir, "day0001.1234.hdf5"))
        sleep(0.1)
        touch(joinpath(output_dir, "day0002.5678.hdf5"))
        restart_file = detect_restart_file(
            ActiveLinkStyle(),
            base_output_dir,
            sort_func = x -> sort(
                x,
                by = f -> stat(joinpath(output_dir, f)).mtime,
                rev = true,
            ),
        )
        @test restart_file == joinpath(output_dir, "day0001.1234.hdf5")
    end
end

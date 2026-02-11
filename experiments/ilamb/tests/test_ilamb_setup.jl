include(joinpath(@__DIR__, "..", "ilamb_setup.jl"))

using Test

"""
    make_dummy_ilamb_diagnostics(dirname)

Make a directory in `tempdir` with the name `dirname` and create two
NetCDF files named "evspsbl_Lmon_CliMA_historical_r1i1p1f1_gn_200003-201902.nc"
and "gpp_Lmon_CliMA_historical_r1i1p1f1_gn_200003-201902.nc" in the directory.

Note that `ilamb_diag_dir` is intentionally not an absolute path.
"""
function make_dummy_ilamb_diagnostics(dirname)
    ilamb_diag_dir = mkdir(dirname)
    evspsbl_filepath = touch(
        joinpath(
            ilamb_diag_dir,
            "evspsbl_Lmon_CliMA_historical_r1i1p1f1_gn_200003-201902.nc",
        ),
    )
    gpp_filepath = touch(
        joinpath(
            ilamb_diag_dir,
            "gpp_Lmon_CliMA_historical_r1i1p1f1_gn_200003-201902.nc",
        ),
    )
    return ilamb_diag_dir, abspath(evspsbl_filepath), abspath(gpp_filepath)
end

function populate_build_dir(build_dir, model_dir)
    isdir(build_dir) || mkdir(build_dir)

    model_name = basename(model_dir)
    touch(joinpath(build_dir, "benchmark_$(model_name).nc"))
    touch(joinpath(build_dir, "benchmark_$(model_name).png"))
    return nothing
end

@testset "Perfect setup for ILAMB using dummy directories and NetCDF files" begin
    tempdir = mktempdir(cleanup = false)
    cd(tempdir) do
        @info "Perfect setup for ILAMB" tempdir

        # Make dummy ILAMB diagnostics
        ilamb_diag_dir, evspsbl_filepath, gpp_filepath =
            make_dummy_ilamb_diagnostics("ilamb_diagnostics")

        # Make dummy ILAMB directory
        ilamb_root_dir = mkdir(joinpath(tempdir, "ilamb"))

        commit_id = "4242424"
        buildkite_id = "1"

        ilamb_build_dir = joinpath(ilamb_root_dir, "_build")
        num_models_to_keep = 2
        setup_run_for_ilamb(
            ilamb_diag_dir,
            ilamb_root_dir,
            ilamb_build_dir,
            buildkite_id,
            commit_id,
            num_models_to_keep,
        )

        # Test existence of directories and symlinks
        MODELS_dir = joinpath(ilamb_root_dir, "MODELS")
        model_name = "$(buildkite_id)_$(Dates.Date(Dates.now()))_$commit_id"
        model_dir = joinpath(MODELS_dir, model_name)
        symlink_evspsbl = joinpath(
            model_dir,
            "evspsbl_Lmon_CliMA_historical_r1i1p1f1_gn_200003-201902.nc",
        )
        symlink_gpp = joinpath(
            model_dir,
            "gpp_Lmon_CliMA_historical_r1i1p1f1_gn_200003-201902.nc",
        )
        @test isdir(MODELS_dir)
        @test isdir(model_dir)
        @test islink(symlink_evspsbl)
        @test islink(symlink_gpp)
        @test readlink(symlink_evspsbl) == evspsbl_filepath
        @test readlink(symlink_gpp) == gpp_filepath

        # Remove files from first setup which simulates deleting old simulation runs
        rm(evspsbl_filepath)
        rm(gpp_filepath)

        # After ILAMB ran, populate the build directory with stuff
        populate_build_dir(ilamb_build_dir, model_dir)

        # Make another ILAMB diagnostics (simulating another run) and set up the
        # dummy run for ILAMB again
        ilamb_diag_dir2, evspsbl_filepath2, gpp_filepath2 =
            make_dummy_ilamb_diagnostics("ilamb_diagnostics1")
        commit_id = "4242425"
        buildkite_id = "2"
        @test_logs (:info, r"Invalid symlink found") match_mode = :any setup_run_for_ilamb(
            ilamb_diag_dir2,
            ilamb_root_dir,
            ilamb_build_dir,
            buildkite_id,
            commit_id,
            num_models_to_keep,
        )

        model_name2 = "$(buildkite_id)_$(Dates.Date(Dates.now()))_$commit_id"
        model_dir2 = joinpath(MODELS_dir, model_name2)
        symlink_evspsbl = joinpath(
            model_dir2,
            "evspsbl_Lmon_CliMA_historical_r1i1p1f1_gn_200003-201902.nc",
        )
        symlink_gpp = joinpath(
            model_dir2,
            "gpp_Lmon_CliMA_historical_r1i1p1f1_gn_200003-201902.nc",
        )
        @test islink(symlink_evspsbl)
        @test islink(symlink_gpp)
        @test readlink(symlink_evspsbl) == evspsbl_filepath2
        @test readlink(symlink_gpp) == gpp_filepath2

        populate_build_dir(ilamb_build_dir, model_dir2)

        # Cleaning the build directory should remove the old png file from
        # populate_build_dir
        files_in_build_dir =
            basename.([
                joinpath(root, file) for
                (root, _, files) in walkdir(ilamb_build_dir) for file in files
            ])
        png_files_in_build_dir =
            filter(filename -> endswith(filename, ".png"), files_in_build_dir)

        @test length(png_files_in_build_dir) == 1
        @test first(png_files_in_build_dir) == "benchmark_$(model_name2).png"

        # Make another ILAMB diagnostics (simulating another run) to test if
        # the first model is deleted
        ilamb_diag_dir3, evspsbl_filepath3, gpp_filepath3 =
            make_dummy_ilamb_diagnostics("ilamb_diagnostics2")
        commit_id = "4242426"
        buildkite_id = "3"
        setup_run_for_ilamb(
            ilamb_diag_dir3,
            ilamb_root_dir,
            ilamb_build_dir,
            buildkite_id,
            commit_id,
            num_models_to_keep,
        )

        model_name3 = "$(buildkite_id)_$(Dates.Date(Dates.now()))_$commit_id"
        model_dir3 = joinpath(MODELS_dir, model_name3)
        populate_build_dir(ilamb_build_dir, model_dir3)

        @test !isdir(model_dir)
        @test isdir(model_dir2)
        @test isdir(model_dir3)

        files_in_build_dir =
            basename.([
                joinpath(root, file) for
                (root, _, files) in walkdir(ilamb_build_dir) for file in files
            ])
        @test length(files_in_build_dir) == 3
        @test isfile(joinpath(ilamb_build_dir, "benchmark_$model_name2.nc"))
        @test isfile(joinpath(ilamb_build_dir, "benchmark_$model_name3.nc"))
        @test isfile(joinpath(ilamb_build_dir, "benchmark_$model_name3.png"))
    end
end

@testset "Model name from simulation" begin
    model_name1 = "1_$(Dates.Date(2010, 3, 2))_a9861b0"
    model_name2 = "1234_$(Dates.Date(2010, 12, 12))_6e04a70"
    model_name3 = "021_$(Dates.Date(1234, 2, 3))_6eac402"
    model_name4 = "CMIP_ENSEMBLE_MODEL"
    model_name5 = "1_2010_1_4_aaaaaaa"
    @test is_model_from_simulation(model_name1)
    @test is_model_from_simulation(model_name2)
    @test is_model_from_simulation(model_name3)
    @test !is_model_from_simulation(model_name4)
    @test !is_model_from_simulation(model_name5)
end

@testset "Delete old runs" begin
    N = 5
    for num_directories in [10, 5, 3]
        temp_dir = mktempdir(cleanup = false)
        @info "Deleting old runs with $num_directories directories" temp_dir
        MODELS_dir = joinpath(temp_dir, "MODELS")
        mkdir(MODELS_dir)

        for i in 1:num_directories
            # Needed for the modified time to not be exactly the same
            model_name = "$(i)_2010-01-01_1234567"
            mkdir(joinpath(MODELS_dir, model_name))
        end

        delete_old_runs(MODELS_dir, N)

        available_dirs = readdir(MODELS_dir, join = true)
        remaining_directories = min(num_directories, N)
        @test length(available_dirs) == remaining_directories
        @test Set(basename.(available_dirs)) == Set([
            "$(i)_2010-01-01_1234567" for
            i in last(1:num_directories, remaining_directories)
        ])
    end

    # If there are extra directories
    temp_dir = mktempdir(cleanup = false)
    @info "Deleting old runs from 4 model directories and 2 directories that should not be removed" temp_dir
    MODELS_dir = joinpath(temp_dir, "MODELS")
    mkdir(MODELS_dir)

    for i in 1:4
        model_name = "$(i)_2010-01-01_1234567"
        mkdir(joinpath(MODELS_dir, model_name))
    end

    mkdir(joinpath(MODELS_dir, "Hello"))
    mkdir(joinpath(MODELS_dir, "World"))

    delete_old_runs(MODELS_dir, N)

    available_dirs = readdir(MODELS_dir, join = true)
    @test length(available_dirs) == 5
    @test Set(basename.(available_dirs)) == union(
        Set(["$(i)_2010-01-01_1234567" for i in 2:4]),
        Set(["Hello", "World"]),
    )
end

@testset "Clean build directory" begin
    # Build directory does not exist
    temp_dir = mktempdir(cleanup = false)
    @info "Cleaning build directory" temp_dir

    MODELS_dir = joinpath(temp_dir, "MODELS")
    mkdir(MODELS_dir)
    for i in 1:3
        model_name = "$(i)_2010-01-01_1234567"
        mkdir(joinpath(MODELS_dir, model_name))
    end

    build_dir = joinpath(temp_dir, "_build")

    @test isnothing(clean_build_dir(MODELS_dir, build_dir))
    @test !isdir(build_dir)


    # Set up fake directory structure
    mkdir(build_dir)
    mkdir(joinpath(build_dir, "1"))
    mkpath(joinpath(build_dir, "2", "3"))
    touch(joinpath(build_dir, "1", "1_2010-01-01_1234567.nc"))
    touch(joinpath(build_dir, "2", "3", "prefix_2_2010-01-01_1234567.nc"))
    touch(joinpath(build_dir, "3_2010-01-01_1234567_suffix.nc"))
    touch(joinpath(build_dir, "misc"))
    touch(joinpath(build_dir, "1", "delete_4_2010-01-01_1234567.nc"))

    clean_build_dir(MODELS_dir, build_dir)

    @test isdir(build_dir)
    @test isdir(joinpath(build_dir, "1"))
    @test isdir(joinpath(build_dir, "2"))
    @test isdir(joinpath(build_dir, "2", "3"))
    @test isfile(joinpath(build_dir, "1", "1_2010-01-01_1234567.nc"))
    @test isfile(
        joinpath(build_dir, "2", "3", "prefix_2_2010-01-01_1234567.nc"),
    )
    @test isfile(build_dir, "3_2010-01-01_1234567_suffix.nc")

    files_in_build_dir =
        basename.([
            joinpath(root, file) for (root, _, files) in walkdir(build_dir) for
            file in files
        ])

    @test length(files_in_build_dir) == 3
    @test "misc" ∉ files_in_build_dir
    @test "delete_4_2010-01-01_1234567.nc" ∉ files_in_build_dir
end

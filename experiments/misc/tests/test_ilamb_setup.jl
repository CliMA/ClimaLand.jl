include(joinpath(@__DIR__, "..", "ilamb_setup.jl"))

using Test

function make_dummy_ilamb_diagnostics(tempdir, dirname)
    ilamb_diag_dir = mkdir(joinpath(tempdir, dirname))
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
    return ilamb_diag_dir, evspsbl_filepath, gpp_filepath
end

@testset "Perfect setup for ILAMB using dummy directories and NetCDF files" begin
    tempdir = mktempdir(cleanup = false)
    @info tempdir

    # Make dummy ILAMB diagnostics
    ilamb_diag_dir, evspsbl_filepath, gpp_filepath =
        make_dummy_ilamb_diagnostics(tempdir, "ilamb_diagnostics")

    # Make dummy ILAMB directory
    ilamb_root_dir = mkdir(joinpath(tempdir, "ilamb"))
    mkdir(joinpath(ilamb_root_dir, "DATA"))

    commit_id = "4242424"
    setup_run_for_ilamb(ilamb_diag_dir, ilamb_root_dir, commit_id)

    MODELS_dir = joinpath(ilamb_root_dir, "MODELS")
    model_dir = joinpath(MODELS_dir, "$(Dates.Date(now()))_$commit_id")
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

    # Remove files from first setup which simulates deleting simulation run if
    # there is not enough space on the cluster
    rm(evspsbl_filepath)
    rm(gpp_filepath)

    # Make another ILAMB diagnostics (simulating another run) and set up the
    # dummy run for ILAMB again
    ilamb_diag_dir2, evspsbl_filepath2, gpp_filepath2 =
        make_dummy_ilamb_diagnostics(tempdir, "ilamb_diagnostics1")
    commit_id = "4242425"
    @test_logs (:info, r"Invalid symlink found") match_mode = :any setup_run_for_ilamb(
        ilamb_diag_dir2,
        ilamb_root_dir,
        commit_id,
    )

    model_dir2 = joinpath(MODELS_dir, "$(Dates.Date(now()))_$commit_id")
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
end

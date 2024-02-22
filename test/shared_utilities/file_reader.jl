#=
    Unit tests for ClimaLand FileReader module
=#

using ClimaLand
using ClimaLand: Bucket, Regridder, FileReader
using ClimaCore
using ClimaCore: Fields, Meshes, Domains, Topologies, Spaces
using ClimaComms
using Test
using Dates
using NCDatasets

# Include testing artifacts from Bucket
albedo_temporal_data = Bucket.cesm2_albedo_dataset_path
albedo_bareground_data = Bucket.bareground_albedo_dataset_path

# Include testing artifacts from test directory
include(joinpath(pkgdir(ClimaLand), "test", "artifacts", "artifacts.jl"))
era5_t2m_sp_u10n_data = era5_t2m_sp_u10n_dataset_path
era5_t2m_sp_u10n_static_data = era5_t2m_sp_u10n_static_dataset_path

comms_ctx = ClimaComms.SingletonCommsContext()

# Use two separate regrid dirs to avoid duplicate filenames since files have same varname
regrid_dir_static = joinpath(pkgdir(ClimaLand), "test", "regridding")
regrid_dir_temporal =
    joinpath(pkgdir(ClimaLand), "test", "regridding", "temporal")

"""
Set NaN and Missing to zero (GPU compatible)
"""
function replace_nan_missing!(field::Fields.Field)
    # For GPU runs, we perform the substitution on the CPU and move back to GPU
    parent_without_NaN_missing =
        replace(Array(parent(field)), NaN => 0, missing => 0)
    ArrayType = ClimaComms.array_type(ClimaComms.device())
    parent(field) .= ArrayType(parent_without_NaN_missing)
end

FT = Float32
@testset "test interpol, FT = $FT" begin
    # Setup
    x1 = FT(0)
    x2 = FT(10)

    f1 = FT(-1)
    f2 = FT(1)

    # Case 1: x1 < x < x2
    x = FT(5)
    @test FileReader.interpol(f1, f2, x - x1, x2 - x1) == 0

    # Case 2: x1 = x < x2
    x = FT(0)
    @test FileReader.interpol(f1, f2, x - x1, x2 - x1) == f1

    # Case 3: x1 < x = x2
    x = FT(10)
    @test FileReader.interpol(f1, f2, x - x1, x2 - x1) == f2
end

@testset "test to_datetime, FT = $FT" begin
    year = 2000
    dt_noleap = DateTimeNoLeap(year)
    @test FileReader.to_datetime(dt_noleap) == DateTime(year)
end

# Add tests which use TempestRemap here -
# TempestRemap is not built on Windows because of NetCDF support limitations
# `PrescribedData` constructors use TR via a call to `hdwrite_regridfile_rll_to_cgll`
if !Sys.iswindows()
    @testset "test PrescribedDataStatic construction, FT = $FT" begin
        # setup for test
        comms_ctx = ClimaComms.SingletonCommsContext()
        radius = FT(6731e3)
        Nq = 4
        domain = ClimaCore.Domains.SphereDomain(radius)
        mesh = Meshes.EquiangularCubedSphere(domain, 4)
        topology = Topologies.Topology2D(comms_ctx, mesh)
        quad = Spaces.Quadratures.GLL{Nq}()
        surface_space_t = Spaces.SpectralElementSpace2D(topology, quad)

        # Loop over files with one and multiple variables
        for (get_infile, varnames) in [
            (albedo_bareground_data, ["sw_alb"]),
            (era5_t2m_sp_u10n_static_data, ["t2m", "sp", "u10n"]),
        ]
            # reset regrid directory between files
            isdir(regrid_dir_static) ? nothing : mkpath(regrid_dir_static)

            ps_data_spatial = FileReader.PrescribedDataStatic{FT}(
                get_infile,
                regrid_dir_static,
                varnames,
                surface_space_t,
            )

            # test that created object exists
            @test @isdefined(ps_data_spatial)
            @test ps_data_spatial isa FileReader.AbstractPrescribedData

            # test fields that we've passed into constructor as args
            @test ps_data_spatial.file_info.infile_path == get_infile()
            @test ps_data_spatial.file_info.regrid_dirpath == regrid_dir_static
            @test ps_data_spatial.file_info.varnames == varnames

            # test fields which are set internally by constructor
            @test ps_data_spatial.file_info.outfile_root == "static_data_cgll"
            @test ps_data_spatial.file_info.all_dates == []
            @test ps_data_spatial.file_info.date_idx0 == []
            rm(regrid_dir_static; recursive = true, force = true)
        end
    end

    @testset "test get_data_at_date for PrescribedDataStatic, FT = $FT" begin
        isdir(regrid_dir_static) ? nothing : mkpath(regrid_dir_static)
        # test `get_data_at_date` with interpolation (default)
        dummy_dates = [DateTime(1999, 1, 1)]
        date_idx0 = Int[1]

        # construct space and dummy field
        comms_ctx = ClimaComms.SingletonCommsContext()
        radius = FT(6731e3)
        Nq = 4
        domain_sphere = ClimaCore.Domains.SphereDomain(radius)
        mesh = Meshes.EquiangularCubedSphere(domain_sphere, 4)
        topology = Topologies.Topology2D(comms_ctx, mesh)
        quad = Spaces.Quadratures.GLL{Nq}()
        surface_space_t = Spaces.SpectralElementSpace2D(topology, quad)
        data_field = ones(surface_space_t) .* FT(0.5)

        # Add 2 variables to `FileInfo`, but only read/write `var1`
        varnames = ["var1", "var2"]
        outfile_root = "static_data_cgll"

        # write data to dummy HDF5 file, which gets read in by `get_data_at_date`
        field = Regridder.write_to_hdf5(
            regrid_dir_static,
            outfile_root,
            Dates.DateTime(0), # dummy date
            data_field,
            varnames[1],
            comms_ctx,
        )

        # create structs manually since we aren't reading from a data file
        file_info = FileReader.FileInfo(
            "",                 # infile_path
            regrid_dir_static,  # regrid_dirpath
            varnames,           # varnames
            outfile_root,       # outfile_root
            dummy_dates,        # all_dates
            date_idx0,          # date_idx0
        )

        prescribed_data_static =
            FileReader.PrescribedDataStatic{typeof(file_info)}(file_info)

        # Read in dummy data that has been written by `write_to_hdf5`
        field_out = FileReader.get_data_at_date(
            prescribed_data_static,
            surface_space_t,
            varnames[1],
        )
        @test field_out == data_field
        rm(regrid_dir_static; recursive = true, force = true)
    end

    @testset "test PrescribedDataTemporal construction, FT = $FT" begin
        # setup for test
        comms_ctx = ClimaComms.SingletonCommsContext()
        radius = FT(6731e3)
        Nq = 4
        domain = ClimaCore.Domains.SphereDomain(radius)
        mesh = Meshes.EquiangularCubedSphere(domain, 4)
        topology = Topologies.Topology2D(comms_ctx, mesh)
        quad = Spaces.Quadratures.GLL{Nq}()
        surface_space_t = Spaces.SpectralElementSpace2D(topology, quad)

        date_idx0 = Int[1]
        date_ref = DateTime(1800, 1, 1)
        t_start = Float64(0)

        # Loop over files with one and multiple variables
        for (get_infile, varnames) in [
            (albedo_temporal_data, ["sw_alb"]),
            (era5_t2m_sp_u10n_data, ["t2m", "sp", "u10n"]),
        ]
            # reset regrid directory between files
            isdir(regrid_dir_temporal) ? nothing : mkpath(regrid_dir_temporal)

            ps_data_temp = FileReader.PrescribedDataTemporal{FT}(
                regrid_dir_temporal,
                get_infile,
                varnames,
                date_ref,
                t_start,
                surface_space_t,
            )

            # test that created object exists
            @test @isdefined(ps_data_temp)
            @test ps_data_temp isa FileReader.AbstractPrescribedData

            # test fields that we've passed into constructor as args
            @test ps_data_temp.file_info.regrid_dirpath == regrid_dir_temporal
            @test ps_data_temp.file_info.varnames == varnames
            @test ps_data_temp.file_info.date_idx0 == date_idx0
            @test ps_data_temp.sim_info.date_ref == date_ref
            @test ps_data_temp.sim_info.t_start == t_start

            # test fields which are set internally by constructor
            @test ps_data_temp.file_info.outfile_root == "temporal_data_cgll"
            @test !isnothing(ps_data_temp.file_info.all_dates)
            @test sort(collect(keys(ps_data_temp.file_states))) ==
                  sort(varnames)

            # check varnames individually
            for varname in varnames
                @test ps_data_temp.file_states[varname].date_idx == date_idx0
                @test ps_data_temp.file_states[varname].date_idx == date_idx0
                @test !isnothing(ps_data_temp.file_states[varname].data_fields)
                @test !isnothing(
                    ps_data_temp.file_states[varname].segment_length,
                )
            end
            rm(regrid_dir_temporal; recursive = true, force = true)
        end
    end

    @testset "test get_data_at_date for PrescribedDataTemporal, FT = $FT" begin
        # test `get_data_at_date` with interpolation (default)
        dummy_dates =
            Vector(range(DateTime(1999, 1, 1); step = Day(1), length = 100))
        date_idx0 = Int[1]
        date_ref = dummy_dates[Int(date_idx0[1]) + 1]
        t_start = Float64(0)
        date0 = date_ref + Dates.Second(t_start)

        # these values give an `interp_fraction` of 0.5 in `interpol` for ease of testing
        segment_length =
            [Int(2) * ((date_ref - dummy_dates[Int(date_idx0[1])]).value)]

        # construct space and dummy field
        comms_ctx = ClimaComms.SingletonCommsContext()
        radius = FT(6731e3)
        Nq = 4
        domain_sphere = ClimaCore.Domains.SphereDomain(radius)
        mesh = Meshes.EquiangularCubedSphere(domain_sphere, 4)
        topology = Topologies.Topology2D(comms_ctx, mesh)
        quad = Spaces.Quadratures.GLL{Nq}()
        surface_space_t = Spaces.SpectralElementSpace2D(topology, quad)
        data_fields = (zeros(surface_space_t), ones(surface_space_t))

        # Add 2 variables to `FileInfo`, but only read/write `var1`
        varnames = ["var1", "var2"]

        # create structs manually since we aren't reading from a data file
        file_info = FileReader.FileInfo(
            "",             # infile_path
            "",             # regrid_dirpath
            varnames,       # varnames
            "",             # outfile_root
            dummy_dates,    # all_dates
            date_idx0,      # date_idx0
        )
        file_state1 = FileReader.FileState(
            data_fields,        # data_fields
            copy(date_idx0),    # date_idx
            segment_length,     # segment_length
        )
        file_state2 = FileReader.FileState(
            data_fields,        # data_fields
            copy(date_idx0),    # date_idx
            segment_length,     # segment_length
        )
        file_states =
            Dict(varnames[1] => file_state1, varnames[2] => file_state2)
        sim_info = FileReader.SimInfo(
            date_ref,       # date_ref
            t_start,        # t_start
        )

        pd_args = (file_info, file_states, sim_info)
        pd_type_args = (
            typeof(file_info),
            typeof(varnames[1]),
            typeof(file_states[varnames[1]]),
            typeof(sim_info),
        )
        prescribed_data =
            FileReader.PrescribedDataTemporal{pd_type_args...}(pd_args...)

        # Note: we expect this to give a warning "No dates available in file..."
        @test FileReader.get_data_at_date(
            prescribed_data,
            surface_space_t,
            varnames[1],
            date0,
        ) == ones(surface_space_t) .* FT(0.5)
    end

    @testset "test next_date_in_file with PrescribedDataTemporal, FT = $FT" begin
        dummy_dates =
            Vector(range(DateTime(1999, 1, 1); step = Day(1), length = 10))
        date_ref = dummy_dates[1]
        t_start = Float64(0)
        date0 = date_ref + Dates.Second(t_start)
        date_idx0 =
            [argmin(abs.(Dates.value(date0) .- Dates.value.(dummy_dates[:])))]
        varname = "var"

        # manually create structs to avoid creating space for outer constructor
        file_info = FileReader.FileInfo(
            "",             # infile_path
            "",             # regrid_dirpath
            [varname],           # varnames
            "",             # outfile_root
            dummy_dates,    # all_dates
            date_idx0,      # date_idx0
        )
        file_state = FileReader.FileState(
            nothing,            # data_fields
            copy(date_idx0),    # date_idx
            Int[],              # segment_length
        )
        file_states = Dict(varname => file_state)
        sim_info = nothing

        pd_args = (file_info, file_states, sim_info)
        pd_type_args = (
            typeof(file_info),
            typeof(varname),
            typeof(file_states[varname]),
            typeof(sim_info),
        )
        prescribed_data =
            FileReader.PrescribedDataTemporal{pd_type_args...}(pd_args...)

        # read in next dates, manually compare to `dummy_dates` array
        idx = date_idx0[1]
        next_date = date0
        for i in 1:(length(dummy_dates) - 1)
            current_date = next_date
            next_date = FileReader.next_date_in_file(prescribed_data)

            @test next_date == dummy_dates[idx + 1]

            prescribed_data.file_states[varname].date_idx[1] += 1
            idx = date_idx0[1] + i
        end
    end

    @testset "test read_data_fields! for single variable, FT = $FT" begin
        # Create regridding directory
        isdir(regrid_dir_temporal) ? nothing : mkpath(regrid_dir_temporal)

        get_infile = albedo_temporal_data
        infile_path = get_infile()
        varname = "sw_alb"
        varnames = [varname]

        # Start with first date in data file
        date0 = NCDataset(infile_path) do ds
            ds["time"][1]
        end
        date0 = FileReader.to_datetime(date0)
        dates = collect(date0:Day(10):(date0 + Day(100))) # includes both endpoints
        date_ref = date0
        t_start = Float64(0)

        comms_ctx = ClimaComms.SingletonCommsContext()
        radius = FT(6731e3)
        Nq = 4
        domain = ClimaCore.Domains.SphereDomain(radius)
        mesh = Meshes.EquiangularCubedSphere(domain, 4)
        topology = Topologies.Topology2D(comms_ctx, mesh)
        quad = Spaces.Quadratures.GLL{Nq}()
        surface_space_t = Spaces.SpectralElementSpace2D(topology, quad)

        prescribed_data = FileReader.PrescribedDataTemporal{FT}(
            regrid_dir_temporal,
            get_infile,
            varnames,
            date_ref,
            t_start,
            surface_space_t,
        )

        # Test each case (1-4)
        current_fields =
            Fields.zeros(FT, surface_space_t), Fields.zeros(FT, surface_space_t)

        # Use this function to reset values between test cases
        function reset_ps_data(ps_data)
            ps_data.file_states[varname].data_fields[1] .= current_fields[1]
            ps_data.file_states[varname].data_fields[2] .= current_fields[2]
            ps_data.file_states[varname].segment_length[1] =
                ps_data.file_states[varname].date_idx[1] =
                    ps_data.file_info.date_idx0[1] = Int(1)
        end

        # Case 1: test date before first date of file, and date aligned with first date of file
        for date in [date0 - Day(1), date0]
            reset_ps_data(prescribed_data)
            prescribed_data_copy = prescribed_data
            FileReader.read_data_fields!(prescribed_data, date, surface_space_t)

            # Test that both data fields store the same data
            # Remove NaNs and missings before comparison
            (; data_fields) = prescribed_data.file_states[varname]

            foreach(replace_nan_missing!, data_fields)

            @test prescribed_data.file_states[varname].data_fields[1] ==
                  prescribed_data.file_states[varname].data_fields[2]
            # date_idx and date_idx0 should be unchanged
            @test prescribed_data.file_states[varname].segment_length == Int[0]
            @test prescribed_data.file_states[varname].date_idx[1] ==
                  prescribed_data_copy.file_states[varname].date_idx[1]
            @test prescribed_data.file_info.date_idx0[1] ==
                  prescribed_data_copy.file_info.date_idx0[1]
        end

        # Case 2: test date after dates in file
        # This case should print 1 warning "this time period is after input data..."
        reset_ps_data(prescribed_data)
        date_after = prescribed_data.file_info.all_dates[end] + Dates.Day(1)
        prescribed_data_copy = prescribed_data
        FileReader.read_data_fields!(
            prescribed_data,
            date_after,
            surface_space_t,
        )

        # Test that both data fields store the same data
        @test prescribed_data.file_states[varname].segment_length == Int[0]

        # Remove NaNs and missings before comparison
        (; data_fields) = prescribed_data.file_states[varname]
        foreach(replace_nan_missing!, data_fields)

        @test prescribed_data.file_states[varname].data_fields[1] ==
              prescribed_data.file_states[varname].data_fields[2]

        # Case 3: loop over simulation dates (general case)
        # This case should print 3 warnings "updating data files: ..."
        data_saved = []
        updating_dates = []

        # Read in data fields over dummy times
        # With the current test setup, read_data_fields! should get called 3 times
        for date in dates
            callback_date = FileReader.next_date_in_file(prescribed_data)

            if (date >= callback_date)
                FileReader.read_data_fields!(
                    prescribed_data,
                    date,
                    surface_space_t,
                )
                push!(
                    data_saved,
                    deepcopy(
                        prescribed_data.file_states[varname].data_fields[1],
                    ),
                )
                push!(updating_dates, deepcopy(date))
            end
        end

        # Replace NaNs and missings for testing
        (; data_fields) = prescribed_data.file_states[varname]
        foreach(replace_nan_missing!, data_saved)

        # Manually read in data from HDF5
        f = prescribed_data.file_states[varname].data_fields[1]
        data_manual = [similar(f), similar(f), similar(f)]
        (; regrid_dirpath, outfile_root, varnames, all_dates) =
            prescribed_data.file_info
        for i in eachindex(data_saved)
            data_manual[i] = Regridder.read_from_hdf5(
                regrid_dirpath,
                outfile_root,
                all_dates[i + 1],
                varnames[1],
                surface_space_t,
            )
            # Replace NaNs and missings for testing comparison
            replace_nan_missing!(data_manual[i])

            @test parent(data_saved[i]) == parent(data_manual[i])
        end

        # Test that the data_saved array was modified
        @test length(data_saved) > 0
        # Check that the final saved date is as expected (midmonth day of last month)
        midmonth_end =
            DateTime(year(dates[end]), month(dates[end]), 15, hour(dates[end]))
        @test updating_dates[end] == midmonth_end

        # Case 4: everything else
        reset_ps_data(prescribed_data)
        prescribed_data.file_states[varname].date_idx[1] =
            prescribed_data.file_info.date_idx0[1] + Int(1)
        date =
            prescribed_data.file_info.all_dates[prescribed_data.file_states[varname].date_idx[1]] -
            Dates.Day(1)

        @test_throws ErrorException FileReader.read_data_fields!(
            prescribed_data,
            date,
            surface_space_t,
        )
        # Delete regridding directory and files
        rm(regrid_dir_temporal; recursive = true, force = true)
    end

    @testset "test read_data_fields! for multiple variables, FT = $FT" begin
        # Create regridding directory
        isdir(regrid_dir_temporal) ? nothing : mkpath(regrid_dir_temporal)

        get_infile = era5_t2m_sp_u10n_data
        infile_path = get_infile()
        varnames = ["t2m", "sp", "u10n"]

        # Start with first date in data file
        date0 = NCDataset(infile_path) do ds
            ds["time"][1]
        end
        date0 = FileReader.to_datetime(date0)
        # Gives dates [00:45, 01:30, 02:15, 03:00, 03:45, 04:30, 05:15, 06:00, 06:45, 07:30]
        dates = collect(date0:Minute(45):(date0 + Minute(450))) # includes both endpoints
        date_ref = date0
        t_start = Float64(0)

        comms_ctx = ClimaComms.SingletonCommsContext()
        radius = FT(6731e3)
        Nq = 4
        domain = ClimaCore.Domains.SphereDomain(radius)
        mesh = Meshes.EquiangularCubedSphere(domain, 4)
        topology = Topologies.Topology2D(comms_ctx, mesh)
        quad = Spaces.Quadratures.GLL{Nq}()
        surface_space_t = Spaces.SpectralElementSpace2D(topology, quad)

        prescribed_data = FileReader.PrescribedDataTemporal{FT}(
            regrid_dir_temporal,
            get_infile,
            varnames,
            date_ref,
            t_start,
            surface_space_t,
        )

        # Test each case (1-4)
        current_fields =
            Fields.zeros(FT, surface_space_t), Fields.zeros(FT, surface_space_t)

        # Use this function to reset values between test cases
        function reset_ps_data(ps_data)
            for varname in varnames
                ps_data.file_states[varname].data_fields[1] .= current_fields[1]
                ps_data.file_states[varname].data_fields[2] .= current_fields[2]
                ps_data.file_states[varname].segment_length[1] =
                    ps_data.file_states[varname].date_idx[1] =
                        ps_data.file_info.date_idx0[1] = Int(1)
            end
        end

        # Case 1: test date before first date of file, and date aligned with first date of file
        for date in [date0 - Day(1), date0]
            reset_ps_data(prescribed_data)
            prescribed_data_copy = prescribed_data
            FileReader.read_data_fields!(prescribed_data, date, surface_space_t)

            # Test that both data fields store the same data
            # Remove NaNs and missings before comparison
            for varname in varnames
                (; data_fields) = prescribed_data.file_states[varname]

                foreach(replace_nan_missing!, data_fields)

                @test prescribed_data.file_states[varname].data_fields[1] ==
                      prescribed_data.file_states[varname].data_fields[2]
                # date_idx and date_idx0 should be unchanged
                @test prescribed_data.file_states[varname].segment_length ==
                      Int[0]
                @test prescribed_data.file_states[varname].date_idx[1] ==
                      prescribed_data_copy.file_states[varname].date_idx[1]
                @test prescribed_data.file_info.date_idx0[1] ==
                      prescribed_data_copy.file_info.date_idx0[1]
            end
        end

        # Case 2: test date after dates in file
        # This case should print 1 warning "this time period is after input data..."
        reset_ps_data(prescribed_data)
        date_after = prescribed_data.file_info.all_dates[end] + Dates.Day(1)
        prescribed_data_copy = prescribed_data
        FileReader.read_data_fields!(
            prescribed_data,
            date_after,
            surface_space_t,
        )

        for varname in varnames
            # Test that both data fields store the same data
            @test prescribed_data.file_states[varname].segment_length == Int[0]
            @test prescribed_data.file_states[varname].data_fields[1] ==
                  prescribed_data.file_states[varname].data_fields[2]
        end

        # Case 3: loop over simulation dates (general case)
        # This case should print 3 warnings "updating data files: ..."
        reset_ps_data(prescribed_data)

        # Read in data fields over dummy times
        # With the current test setup, read_data_fields! should get called 3 times
        data_saved = Dict{String, Vector{Fields.Field}}(
            varnames[1] => [],
            varnames[2] => [],
            varnames[3] => [],
        )
        updating_dates = []
        for date in dates
            callback_date = FileReader.next_date_in_file(prescribed_data)

            if (date >= callback_date)
                FileReader.read_data_fields!(
                    prescribed_data,
                    date,
                    surface_space_t,
                )
                for varname in varnames
                    push!(
                        data_saved[varname],
                        deepcopy(
                            prescribed_data.file_states[varname].data_fields[1],
                        ),
                    )
                    # Replace NaNs and missings for testing
                    replace_nan_missing!(data_saved[varname][end])
                end

                push!(updating_dates, callback_date)
            end
        end

        # Manually read in data from HDF5
        data_manual = Dict{String, Vector{Fields.Field}}(
            varnames[1] => [],
            varnames[2] => [],
            varnames[3] => [],
        )
        (; regrid_dirpath, outfile_root, varnames, all_dates) =
            prescribed_data.file_info

        for i in eachindex(updating_dates)
            for varname in varnames
                push!(
                    data_manual[varname],
                    Regridder.read_from_hdf5(
                        regrid_dirpath,
                        outfile_root,
                        updating_dates[i],
                        varname,
                        surface_space_t,
                    ),
                )
                # Replace NaNs and missings for testing comparison
                replace_nan_missing!(data_manual[varname][end])
            end
        end

        # Test read_data_fields and manually read data are the same
        for varname in varnames
            @test all(
                parent(data_saved[varname]) .== parent(data_manual[varname]),
            )
        end

        # Test that the data_saved array was modified
        @test length(data_saved) > 0
        # Check that the final saved date is as expected (midmonth day of last month)
        midmonth_end = DateTime(
            year(dates[end]),
            month(dates[end]),
            day(dates[end]),
            hour(dates[end]),
        )
        @test updating_dates[end] == midmonth_end

        # Case 4: everything else
        reset_ps_data(prescribed_data)
        for varname in varnames
            prescribed_data.file_states[varname].date_idx[1] =
                prescribed_data.file_info.date_idx0[1] + Int(1)
            date =
                prescribed_data.file_info.all_dates[prescribed_data.file_states[varname].date_idx[1]] -
                Dates.Day(1)

            @test_throws ErrorException FileReader.read_data_fields!(
                prescribed_data,
                date,
                surface_space_t,
            )
        end
        # Delete regridding directory and files
        rm(regrid_dir_temporal; recursive = true, force = true)
    end
end

# Delete testing directory and files
rm(regrid_dir_static; recursive = true, force = true)
rm(regrid_dir_temporal; recursive = true, force = true)

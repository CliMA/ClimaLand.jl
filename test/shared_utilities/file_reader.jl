#=
    Unit tests for ClimaLSM FileReader module
=#

using ClimaLSM
using ClimaLSM: Bucket, Regridder, FileReader
using ClimaCore
using ClimaCore: Fields, Meshes, Domains, Topologies, Spaces
using ClimaComms
using Test
using Dates
using NCDatasets

albedo_temporal_data = Bucket.cesm2_albedo_dataset_path()
albedo_bareground_data = Bucket.bareground_albedo_dataset_path()

comms_ctx = ClimaComms.SingletonCommsContext()

# Use two separate regrid dirs to avoid duplicate filenames since files have same varname
regrid_dir_static = joinpath(pkgdir(ClimaLSM), "test", "regridding")
regrid_dir_temporal =
    joinpath(pkgdir(ClimaLSM), "test", "regridding", "temporal")
isdir(regrid_dir_static) ? nothing : mkpath(regrid_dir_static)
isdir(regrid_dir_temporal) ? nothing : mkpath(regrid_dir_temporal)


FT = Float32
@testset "test interpol, FT = $FT" begin
    # Setup
    t1 = FT(0)
    t2 = FT(10)

    f1 = FT(-1)
    f2 = FT(1)

    # Case 1: t1 < t < t2
    t = FT(5)
    @test FileReader.interpol(f1, f2, t - t1, t2 - t1) == 0

    # Case 2: t1 = t < t2
    t = FT(0)
    @test FileReader.interpol(f1, f2, t - t1, t2 - t1) == f1

    # Case 3: t1 < t = t2
    t = FT(10)
    @test FileReader.interpol(f1, f2, t - t1, t2 - t1) == f2
end

@testset "test interpolate_data, FT = $FT" begin
    # test interpolate_data with interpolation (default)
    dummy_dates =
        Vector(range(DateTime(1999, 1, 1); step = Day(1), length = 100))
    date_idx0 = Int[1]
    date_ref = dummy_dates[Int(date_idx0[1]) + 1]
    t_start = FT(0)
    date0 = date_ref + Dates.Second(t_start)

    # these values give an `interp_fraction` of 0.5 in `interpol` for ease of testing
    segment_length =
        [Int(2) * ((date_ref - dummy_dates[Int(date_idx0[1])]).value)]

    radius = FT(6731e3)
    Nq = 4
    domain_sphere = ClimaCore.Domains.SphereDomain(radius)
    mesh = Meshes.EquiangularCubedSphere(domain_sphere, 4)
    topology = Topologies.Topology2D(mesh)
    quad = Spaces.Quadratures.GLL{Nq}()
    surface_space_t = Spaces.SpectralElementSpace2D(topology, quad)
    data_fields = (zeros(surface_space_t), ones(surface_space_t))

    # create structs manually since we aren't reading from a data file
    file_info = FileReader.FileInfo(
        "",             # infile_path
        "",             # regrid_dirpath
        "",             # varname
        "",             # outfile_root
        dummy_dates,    # all_dates
        date_idx0,      # date_idx0
    )
    file_state = FileReader.FileState(
        data_fields,        # data_fields
        copy(date_idx0),    # date_idx
        segment_length,     # segment_length
    )
    sim_info = FileReader.SimInfo(
        date_ref,       # date_ref
        t_start,        # t_start
    )

    pd_args = (file_info, file_state, sim_info)
    prescribed_data =
        FileReader.PrescribedDataTemporal{typeof.(pd_args)...}(pd_args...)

    # Note: we expect this to give a warning "No dates available in file..."
    @test FileReader.interpolate_data(
        prescribed_data,
        date0,
        surface_space_t,
    ) == ones(surface_space_t) .* FT(0.5)
end

@testset "test to_datetime, FT = $FT" begin
    year = 2000
    dt_noleap = DateTimeNoLeap(year)
    @test FileReader.to_datetime(dt_noleap) == DateTime(year)
end

@testset "test next_date_in_file, FT = $FT" begin
    dummy_dates =
        Vector(range(DateTime(1999, 1, 1); step = Day(1), length = 10))
    date_ref = dummy_dates[1]
    t_start = FT(0)
    date0 = date_ref + Dates.Second(t_start)
    date_idx0 =
        [argmin(abs.(Dates.value(date0) .- Dates.value.(dummy_dates[:])))]

    # manually create structs to avoid creating space for outer constructor
    file_info = FileReader.FileInfo(
        "",             # infile_path
        "",             # regrid_dirpath
        "",             # varname
        "",             # outfile_root
        dummy_dates,    # all_dates
        date_idx0,      # date_idx0
    )
    file_state = FileReader.FileState(
        nothing,            # data_fields
        copy(date_idx0),    # date_idx
        Int[],              # segment_length
    )
    sim_info = nothing

    pd_args = (file_info, file_state, sim_info)
    prescribed_data =
        FileReader.PrescribedDataTemporal{typeof.(pd_args)...}(pd_args...)

    # read in next dates, manually compare to `dummy_dates` array
    idx = date_idx0[1]
    next_date = date0
    for i in 1:(length(dummy_dates) - 1)
        current_date = next_date
        next_date = FileReader.next_date_in_file(prescribed_data)

        @test next_date == dummy_dates[idx + 1]

        prescribed_data.file_state.date_idx[1] += 1
        idx = date_idx0[1] + i
    end
end

@testset "test PrescribedDataStatic construction, FT = $FT" begin
    infile_path = albedo_temporal_data
    varname = "sw_alb"
    ps_data_spatial =
        FileReader.PrescribedDataStatic(infile_path, regrid_dir_static, varname)

    # test that created object exists
    @test @isdefined(ps_data_spatial)
    @test ps_data_spatial isa FileReader.AbstractPrescribedData

    # test fields that we've passed into constructor as args
    @test ps_data_spatial.file_info.infile_path == infile_path
    @test ps_data_spatial.file_info.regrid_dirpath == regrid_dir_static
    @test ps_data_spatial.file_info.varname == varname

    # test fields which are set internally by constructor
    @test ps_data_spatial.file_info.outfile_root == ""
    @test ps_data_spatial.file_info.all_dates == []
    @test ps_data_spatial.file_info.date_idx0 == []
end

# Add tests which use TempestRemap here -
# TempestRemap is not built on Windows because of NetCDF support limitations
# `PrescribedDataTemporal` uses TR via a call to `hdwrite_regridfile_rll_to_cgll`
if !Sys.iswindows()
    @testset "test PrescribedDataTemporal construction, FT = $FT" begin
        # setup for test
        radius = FT(6731e3)
        Nq = 4
        domain = ClimaCore.Domains.SphereDomain(radius)
        mesh = Meshes.EquiangularCubedSphere(domain, 4)
        topology = Topologies.Topology2D(mesh)
        quad = Spaces.Quadratures.GLL{Nq}()
        surface_space_t = Spaces.SpectralElementSpace2D(topology, quad)

        infile_path = albedo_temporal_data
        varname = "sw_alb"
        date_idx0 = Int[1]
        date_ref = DateTime(1800, 1, 1)
        t_start = FT(0)

        ps_data_temp = FileReader.PrescribedDataTemporal{FT}(
            regrid_dir_temporal,
            infile_path,
            varname,
            date_ref,
            t_start,
            surface_space_t,
        )

        # test that created object exists
        @test @isdefined(ps_data_temp)
        @test ps_data_temp isa FileReader.AbstractPrescribedData

        # test fields that we've passed into constructor as args
        @test ps_data_temp.file_info.regrid_dirpath == regrid_dir_temporal
        @test ps_data_temp.file_info.varname == varname
        @test ps_data_temp.file_info.date_idx0 == date_idx0
        @test ps_data_temp.file_state.date_idx == date_idx0
        @test ps_data_temp.sim_info.date_ref == date_ref
        @test ps_data_temp.sim_info.t_start == t_start

        # test fields which are set internally by constructor
        @test ps_data_temp.file_info.outfile_root != ""
        @test !isnothing(ps_data_temp.file_info.all_dates)
        @test !isnothing(ps_data_temp.file_state.data_fields)
        @test !isnothing(ps_data_temp.file_state.segment_length)
    end

    @testset "test read_data_fields!, FT = $FT" begin
        # setup for test
        datafile_rll = albedo_temporal_data
        varname = "sw_alb"

        # Start with first date in data file
        date0 = NCDataset(datafile_rll) do ds
            ds["time"][1]
        end
        date0 = FileReader.to_datetime(date0)
        dates = collect(date0:Day(10):(date0 + Day(100))) # includes both endpoints
        date_ref = date0
        t_start = FT(0)

        radius = FT(6731e3)
        Nq = 4
        domain = ClimaCore.Domains.SphereDomain(radius)
        mesh = Meshes.EquiangularCubedSphere(domain, 4)
        topology = Topologies.DistributedTopology2D(
            comms_ctx,
            mesh,
            Topologies.spacefillingcurve(mesh),
        )
        quad = Spaces.Quadratures.GLL{Nq}()
        surface_space_t = Spaces.SpectralElementSpace2D(topology, quad)

        prescribed_data = FileReader.PrescribedDataTemporal{FT}(
            regrid_dir_temporal,
            albedo_temporal_data,
            varname,
            date_ref,
            t_start,
            surface_space_t,
        )

        # Test each case (1-4)
        current_fields =
            Fields.zeros(FT, surface_space_t), Fields.zeros(FT, surface_space_t)

        # Use this function to reset values between test cases
        function reset_ps_data(ps_data)
            ps_data.file_state.data_fields[1] .= current_fields[1]
            ps_data.file_state.data_fields[2] .= current_fields[2]
            ps_data.file_state.segment_length[1] =
                ps_data.file_state.date_idx[1] =
                    ps_data.file_info.date_idx0[1] = Int(1)
        end

        # Case 1: test initial date - aligned with first date of albedo file
        # This case should not print any warnings
        reset_ps_data(prescribed_data)
        prescribed_data_copy = prescribed_data
        FileReader.read_data_fields!(prescribed_data, date0, surface_space_t)

        # Test that both data fields store the same data
        # Remove NaNs and missings before comparison
        (; data_fields) = prescribed_data.file_state
        replace!(parent(data_fields[1]), NaN => 0)
        replace!(parent(data_fields[1]), missing => 0)
        replace!(parent(data_fields[2]), NaN => 0)
        replace!(parent(data_fields[2]), missing => 0)

        @test prescribed_data.file_state.data_fields[1] ==
              prescribed_data.file_state.data_fields[2]
        # date_idx and date_idx0 should be unchanged
        @test prescribed_data.file_state.segment_length == Int[0]
        @test prescribed_data.file_state.date_idx[1] ==
              prescribed_data_copy.file_state.date_idx[1]
        @test prescribed_data.file_info.date_idx0[1] ==
              prescribed_data_copy.file_info.date_idx0[1]


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
        @test prescribed_data.file_state.segment_length == Int[0]

        # Remove NaNs and missings before comparison
        replace!(parent(prescribed_data.file_state.data_fields[1]), NaN => 0)
        replace!(
            parent(prescribed_data.file_state.data_fields[1]),
            missing => 0,
        )
        replace!(parent(prescribed_data.file_state.data_fields[2]), NaN => 0)
        replace!(
            parent(prescribed_data.file_state.data_fields[2]),
            missing => 0,
        )
        @test prescribed_data.file_state.data_fields[1] ==
              prescribed_data.file_state.data_fields[2]

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
                    deepcopy(prescribed_data.file_state.data_fields[1]),
                )
                push!(updating_dates, deepcopy(date))
            end
        end

        # Replace NaNs and missings for testing
        for ds in data_saved
            replace!(parent(ds), NaN => 0)
            replace!(parent(ds), missing => 0)
        end

        # Manually read in data from HDF5
        f = prescribed_data.file_state.data_fields[1]
        data_manual = [similar(f), similar(f), similar(f)]
        (; regrid_dirpath, outfile_root, varname, all_dates) =
            prescribed_data.file_info
        for i in eachindex(data_saved)
            data_manual[i] = Regridder.swap_space!(
                Regridder.read_from_hdf5(
                    regrid_dirpath,
                    outfile_root,
                    all_dates[i + 1],
                    varname,
                    comms_ctx,
                ),
                surface_space_t,
            )
            # Replace NaNs and missings for testing comparison
            replace!(parent(data_manual[i]), NaN => 0)
            replace!(parent(data_manual[i]), missing => 0)

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
        prescribed_data.file_state.date_idx[1] =
            prescribed_data.file_info.date_idx0[1] + Int(1)
        date =
            prescribed_data.file_info.all_dates[prescribed_data.file_state.date_idx[1]] -
            Dates.Day(1)

        @test_throws ErrorException FileReader.read_data_fields!(
            prescribed_data,
            date,
            surface_space_t,
        )
    end
end

# Delete testing directory and files
rm(regrid_dir_static; recursive = true, force = true)
rm(regrid_dir_temporal; recursive = true, force = true)

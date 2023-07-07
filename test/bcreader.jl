#=
    Unit tests for ClimaLSM BCReader module
=#

using ClimaLSM
using ClimaLSM: Bucket, Regridder, BCReader
using ClimaCore: Fields, Meshes, Domains, Topologies, Spaces
using ClimaComms
using Test
using Dates
using NCDatasets
using CFTime

include(joinpath(pkgdir(ClimaLSM), "src/Bucket/artifacts/artifacts.jl"))
albedo_temporal_data = cesm2_albedo_dataset_path()
mask_data = mask_dataset_path()

comms_ctx = ClimaComms.SingletonCommsContext()
regrid_dir = joinpath(pkgdir(ClimaLSM), "test", "regridding")
isdir(regrid_dir) ? nothing : mkpath(regrid_dir)


for FT in (Float64, Float32)
    @testset "test to_datetime for FT=$FT" begin
        year = 2000
        dt_noleap = DateTimeNoLeap(year)
        @test BCReader.to_datetime(dt_noleap) == DateTime(year)
    end

    @testset "test datetime_to_strdate for FT=$FT" begin
        date1 = DateTime(1950, 12, 31)
        @test BCReader.datetime_to_strdate(date1) == "19501231"

        date2 = date1 + Dates.Day(1)
        @test BCReader.datetime_to_strdate(date2) == "19510101"
    end

    @testset "test next_date_in_file for FT=$FT" begin
        dummy_dates =
            Vector(range(DateTime(1999, 1, 1); step = Day(1), length = 10))
        date0 = dummy_dates[1]
        segment_idx0 = [
            argmin(
                abs.(
                    parse(FT, BCReader.datetime_to_strdate(date0)) .-
                    parse.(FT, BCReader.datetime_to_strdate.(dummy_dates[:]))
                ),
            ),
        ]

        bcf_info = BCReader.BCFileInfo{FT}(
            "",                                 # bcfile_dir
            comms_ctx,                          # comms_ctx
            "",                                 # hd_outfile_root
            "",                                 # varname
            dummy_dates,                        # all_dates
            nothing,                            # monthly_fields
            nothing,                            # scaling_function
            nothing,                            # land_fraction
            deepcopy(segment_idx0),             # segment_idx
            segment_idx0,                       # segment_idx0
            Int[],                              # segment_length
            false,                              # interpolate_daily
        )

        idx = segment_idx0[1]
        current_date = date0
        next_date = BCReader.next_date_in_file(bcf_info)
        @test current_date == dummy_dates[idx]

        for i in 1:(length(dummy_dates) - 2)
            current_date = next_date
            bcf_info.segment_idx[1] += Int(1)
            next_date = BCReader.next_date_in_file(bcf_info)
            idx = segment_idx0[1] + i
            @test next_date == dummy_dates[idx + 1]
        end
    end

    @testset "test interpolate_midmonth_data for FT=$FT" begin
        # test interpolate_midmonth_data with interpolation
        interpolate_daily = true
        dummy_dates =
            Vector(range(DateTime(1999, 1, 1); step = Day(1), length = 100))
        segment_idx0 = [Int(1)]

        # these values give an `interp_fraction` of 0.5 in `interpol` for ease of testing
        date0 = dummy_dates[Int(segment_idx0[1] + 1)]
        segment_length =
            [Int(2) * ((date0 - dummy_dates[Int(segment_idx0[1])]).value)]

        radius = FT(6731e3)
        Nq = 4
        domain_sphere = Domains.SphereDomain(radius)
        mesh = Meshes.EquiangularCubedSphere(domain_sphere, 4)
        topology = Topologies.Topology2D(mesh)
        quad = Spaces.Quadratures.GLL{Nq}()
        boundary_space_t = Spaces.SpectralElementSpace2D(topology, quad)
        monthly_fields = (zeros(boundary_space_t), ones(boundary_space_t))

        bcf_info_interp = BCReader.BCFileInfo{FT}(
            "",                                 # bcfile_dir
            comms_ctx,                          # comms_ctx
            "",                                 # hd_outfile_root
            "",                                 # varname
            dummy_dates,                        # all_dates
            monthly_fields,                     # monthly_fields
            nothing,                            # scaling_function
            nothing,                            # land_fraction
            deepcopy(segment_idx0),             # segment_idx
            segment_idx0,                       # segment_idx0
            segment_length,                     # segment_length
            interpolate_daily,                  # interpolate_daily
        )
        @test BCReader.interpolate_midmonth_data(date0, bcf_info_interp) ==
              ones(boundary_space_t) .* FT(0.5)

        # test interpolate_midmonth_data without interpolation
        interpolate_daily = false

        bcf_info_no_interp = BCReader.BCFileInfo{FT}(
            "",                                 # bcfile_dir
            comms_ctx,                          # comms_ctx
            "",                                 # hd_outfile_root
            "",                                 # varname
            dummy_dates,                        # all_dates
            monthly_fields,                     # monthly_fields
            nothing,                            # scaling_function
            nothing,                            # land_fraction
            deepcopy(segment_idx0),             # segment_idx
            segment_idx0,                       # segment_idx0
            segment_length,                     # segment_length
            interpolate_daily,                  # interpolate_daily
        )
        @test BCReader.interpolate_midmonth_data(date0, bcf_info_no_interp) ==
              monthly_fields[1]
    end

    # Add tests which use TempestRemap here -
    # TempestRemap is not built on Windows because of NetCDF support limitations
    # `bcf_info_init` uses TR via a call to `hdwrite_regridfile_rll_to_cgll`
    if !Sys.iswindows()
        @testset "test bcf_info_init for FT=$FT" begin
            # setup for test
            radius = FT(6731e3)
            Nq = 4
            domain = Domains.SphereDomain(radius)
            mesh = Meshes.EquiangularCubedSphere(domain, 4)
            topology = Topologies.Topology2D(mesh)
            quad = Spaces.Quadratures.GLL{Nq}()
            boundary_space_t = Spaces.SpectralElementSpace2D(topology, quad)
            land_fraction_t = Fields.zeros(boundary_space_t)

            datafile_rll = mask_data
            varname = "LSMASK"
            mono = true

            bcf_info = BCReader.bcfile_info_init(
                FT,
                regrid_dir,
                datafile_rll,
                varname,
                boundary_space_t,
                comms_ctx,
                segment_idx0 = [Int(1309)],
                land_fraction = land_fraction_t,
                mono = mono,
            )

            # test that created object exists and has correct components
            @test @isdefined(bcf_info)
            @test all(parent(bcf_info.land_fraction) .== 0)
        end

        @testset "test update_midmonth_data! for FT=$FT" begin
            # setup for test
            datafile_rll = albedo_temporal_data
            varname = "sw_alb"

            # Start with first date in data file
            date0 = BCReader.to_datetime(NCDataset(datafile_rll)["time"][1])
            dates = collect(date0:Day(10):(date0 + Day(100))) # includes both endpoints

            radius = FT(6731e3)
            Nq = 4
            domain = Domains.SphereDomain(radius)
            mesh = Meshes.EquiangularCubedSphere(domain, 4)
            topology = Topologies.DistributedTopology2D(
                comms_ctx,
                mesh,
                Topologies.spacefillingcurve(mesh),
            )
            quad = Spaces.Quadratures.GLL{Nq}()
            boundary_space_t = Spaces.SpectralElementSpace2D(topology, quad)

            land_fraction_t = Fields.ones(boundary_space_t)

            bcf_info = BCReader.bcfile_info_init(
                FT,
                regrid_dir,
                datafile_rll,
                varname,
                boundary_space_t,
                comms_ctx,
                interpolate_daily = true,
                land_fraction = land_fraction_t,
            )

            # General case: loop over simulation dates
            data_saved = []
            updating_dates = []

            for date in dates
                callback_date = BCReader.next_date_in_file(bcf_info)

                if (date >= callback_date)
                    BCReader.update_midmonth_data!(date, bcf_info)
                    push!(data_saved, deepcopy(bcf_info.monthly_fields[1]))
                    push!(updating_dates, deepcopy(date))
                end
            end

            # Test that the data_saved array was modified
            @test length(data_saved) > 0
            # Check that the final saved date is as expected (midmonth day of last month)
            midmonth_end = DateTime(
                year(dates[end]),
                month(dates[end]),
                15,
                hour(dates[end]),
            )
            @test updating_dates[end] == midmonth_end

            # Test warning/error cases
            current_fields = Fields.zeros(FT, boundary_space_t),
            Fields.zeros(FT, boundary_space_t)

            # Use this function to reset values between test cases
            function reset_bcf_info(bcf_info)
                bcf_info.monthly_fields[1] .= current_fields[1]
                bcf_info.monthly_fields[2] .= current_fields[2]
                bcf_info.segment_length[1] =
                    bcf_info.segment_idx[1] = bcf_info.segment_idx0[1] = Int(1)
            end

            # Case 1: test initial date - aligned with first date of albedo file
            reset_bcf_info(bcf_info)
            bcf_info_copy = bcf_info
            BCReader.update_midmonth_data!(date0, bcf_info)

            if !any(isnan.(parent(bcf_info.monthly_fields[1])))
                @test bcf_info.monthly_fields[1] == bcf_info.monthly_fields[2]
            end
            # bcf_info.segment_idx shouldn't get modified even if midmonth_idx modified
            @test bcf_info.segment_length == [Int(0)]
            @test bcf_info.segment_idx[1] == bcf_info_copy.segment_idx[1]
            @test bcf_info.segment_idx0[1] == bcf_info_copy.segment_idx0[1]


            # Case 2: test date after dates in file
            reset_bcf_info(bcf_info)
            date_after = bcf_info.all_dates[end] + Dates.Day(1)
            bcf_info_copy = bcf_info
            BCReader.update_midmonth_data!(date_after, bcf_info)

            @test bcf_info.segment_length == [Int(0)]
            if !any(isnan.(parent(bcf_info.monthly_fields[1])))
                @test bcf_info.monthly_fields[1] == bcf_info.monthly_fields[2]
            end

            # Case 3: closer initial indices in BC file matching this date (?)
            reset_bcf_info(bcf_info)
            date = DateTime(
                bcf_info.all_dates[bcf_info.segment_idx[1] + 1] + Dates.Day(3),
            )
            BCReader.update_midmonth_data!(date, bcf_info)

            nearest_idx = argmin(
                abs.(
                    parse(FT, BCReader.datetime_to_strdate(date)) .-
                    parse.(
                        FT,
                        BCReader.datetime_to_strdate.(bcf_info.all_dates[:]),
                    )
                ),
            )
            @test bcf_info.segment_idx[1] == nearest_idx

            # Case 4: everything else
            reset_bcf_info(bcf_info)
            bcf_info.segment_idx[1] = bcf_info.segment_idx0[1] + Int(1)
            date = bcf_info.all_dates[bcf_info.segment_idx[1]] - Dates.Day(1)

            @test_throws ErrorException BCReader.update_midmonth_data!(
                date,
                bcf_info,
            )

        end
    end
end

# delete testing directory and files
rm(regrid_dir; recursive = true, force = true)

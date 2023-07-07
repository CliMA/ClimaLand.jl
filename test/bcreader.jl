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

# TODO add loop over FTs

FT = Float64


"""
these pass :)
@testset "test to_datetime" begin
    year = 2000
    dt_noleap = DateTimeNoLeap(year)
    @test BCReader.to_datetime(dt_noleap) == DateTime(year)
end

@testset "test datetime_to_strdate" begin
    date1 = DateTime(1950, 12, 31)
    @test BCReader.datetime_to_strdate(date1) == "19501231"

    date2 = date1 + Dates.Day(1)
    @test BCReader.datetime_to_strdate(date2) == "19510101"
end

@testset "test next_date_in_file for FT=$FT" begin
    dummy_dates = Vector(range(DateTime(1999, 1, 1); step = Day(1), length = 10))
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

@testset "test interpolate_midmonth_to_daily for FT=$FT" begin
    # test interpolate_midmonth_to_daily with interpolation
    interpolate_daily = true
    dummy_dates = Vector(range(DateTime(1999, 1, 1); step = Day(1), length = 100))
    segment_idx0 = [Int(1)]

    # these values give an `interp_fraction` of 0.5 in `interpol` for ease of testing
    date0 = dummy_dates[Int(segment_idx0[1] + 1)]
    segment_length = [Int(2) * ((date0 - dummy_dates[Int(segment_idx0[1])]).value)]

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
    @test BCReader.interpolate_midmonth_to_daily(date0, bcf_info_interp) == ones(boundary_space_t) .* FT(0.5)

    # test interpolate_midmonth_to_daily without interpolation
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
    @test BCReader.interpolate_midmonth_to_daily(date0, bcf_info_no_interp) == monthly_fields[1]
end
"""

# # Add tests which use TempestRemap here -
# # TempestRemap is not built on Windows because of NetCDF support limitations
# # `bcf_info_init` uses TR via a call to `hdwrite_regridfile_rll_to_cgll`
# if !Sys.iswindows()
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

#     @testset "test update_midmonth_data! for FT=$FT" begin
        # setup for test
        # datafile_rll = albedo_temporal_data
        # varname = "sw_alb"

        # # Start with first date in data file
        # date0 = BCReader.to_datetime(NCDataset(datafile_rll)["time"][1])
        # dates = collect(date0:Day(10):(date0 + Day(100))) # includes both endpoints
        # # date0 = date1 = DateTime(1979, 01, 01, 01, 00, 00)
        # # date = DateTime(1979, 01, 01, 00, 00, 00)
        # # tspan = (Int(1), Int(90 * 86400)) # Jan-Mar
        # # Δt = Int(1 * 3600)

        # radius = FT(6731e3)
        # Nq = 4
        # domain = Domains.SphereDomain(radius)
        # mesh = Meshes.EquiangularCubedSphere(domain, 4)
        # topology = Topologies.DistributedTopology2D(comms_ctx, mesh, Topologies.spacefillingcurve(mesh))
        # quad = Spaces.Quadratures.GLL{Nq}()
        # boundary_space_t = Spaces.SpectralElementSpace2D(topology, quad)

        # land_fraction_t = Fields.ones(boundary_space_t) # TODO test real land mask?
        # # dummy_data = (; test_data = zeros(axes(land_fraction_t)))

        # bcf_info = BCReader.bcfile_info_init(
        #     FT,
        #     regrid_dir,
        #     datafile_rll,
        #     varname,
        #     boundary_space_t,
        #     comms_ctx,
        #     interpolate_daily = true,
        #     land_fraction = land_fraction_t,
        # )

        # # Test initial date - aligned with first date of albedo file
        # bcf_info_copy = bcf_info
        # BCReader.update_midmonth_data!(date0, bcf_info)

        # @test bcf_info.monthly_fields[1] == bcf_info.monthly_fields[2]
        # @test bcf_info.segment_length == Int(0)

        # # TODO not sure if this is right - does bcf_info field get modified when midmonth_idx modified?
        # @test bcf_info.segment_idx[1] == bcf_info_copy.segment_idx[1] - 1
        # @test bcf_info.segment_idx0[1] == bcf_info_copy.segment_idx0[1]


        # Test following dates - not all aligned with dates of albedo file

"""
        # TODO adapt test from coupler to ClimaLSM
        # dates = (; date = [date], date0 = [date0], date1 = [date1])
        # SST_all = []
        # updating_dates = []

        # # step in time
        # walltime = @elapsed for t in ((tspan[1] + Δt):Δt:tspan[end])
        #     cs_t.dates.date[1] = current_date(cs_t, t) # if not global, `date`` is not updated. Check that this still runs when distributed.

        #     model_date = cs_t.dates.date[1]
        #     callback_date = BCReader.next_date_in_file(bcf_info)

        #     # TODO investigate if macro would be faster here
        #     if (model_date >= callback_date)
        #         BCReader.update_midmonth_data!(model_date, bcf_info)
        #         push!(SST_all, deepcopy(bcf_info.monthly_fields[1]))
        #         push!(updating_dates, deepcopy(model_date))
        #     end

        # end

        # # test if the SST field was modified
        # @test SST_all[end] !== SST_all[end - 1]
        # # check that the final file date is as expected
        # @test Date(updating_dates[end]) == Date(1979, 03, 16)

        # # test warning/error cases
        # current_fields = Fields.zeros(FT, boundary_space_t), Fields.zeros(FT, boundary_space_t)

        # # use this function to reset values between test cases
        # function reset_bcf_info(bcf_info)
        #     bcf_info.monthly_fields[1] .= current_fields[1]
        #     bcf_info.monthly_fields[2] .= current_fields[2]
        #     bcf_info.segment_length[1] = Int(1)
        # end

        # hd_outfile_root = varname * "_cgll"

        # #  case 1: segment_idx == segment_idx0, date < all_dates[segment_idx]
        # bcf_info.segment_idx[1] = bcf_info.segment_idx0[1]
        # date = DateTime(bcf_info.all_dates[bcf_info.segment_idx[1]] - Dates.Day(1))
        # BCReader.update_midmonth_data!(date, bcf_info)

        # @test bcf_info.monthly_fields[1] == bcf_info.scaling_function(
        #     Regridder.read_from_hdf5(
        #         regrid_dir,
        #         hd_outfile_root,
        #         bcf_info.all_dates[Int(bcf_info.segment_idx0[1])],
        #         varname,
        #         comms_ctx,
        #     ),
        #     bcf_info,
        # )
        # @test bcf_info.monthly_fields[2] == bcf_info.monthly_fields[1]
        # @test bcf_info.segment_length[1] == Int(0)

        # #  case 2: date > all_dates[end - 1]
        # reset_bcf_info(bcf_info)
        # date = DateTime(bcf_info.all_dates[end - 1] + Dates.Day(1))
        # BCReader.update_midmonth_data!(date, bcf_info)

        # @test bcf_info.monthly_fields[1] == bcf_info.scaling_function(
        #     Regridder.read_from_hdf5(
        #         regrid_dir,
        #         hd_outfile_root,
        #         bcf_info.all_dates[Int(length(bcf_info.all_dates))],
        #         varname,
        #         comms_ctx,
        #     ),
        #     bcf_info,
        # )
        # @test bcf_info.monthly_fields[2] == bcf_info.monthly_fields[1]
        # @test bcf_info.segment_length[1] == Int(0)

        # #  case 3: date - all_dates[segment_idx + 1] > 2
        # reset_bcf_info(bcf_info)
        # date = DateTime(bcf_info.all_dates[bcf_info.segment_idx[1] + 1] + Dates.Day(3))
        # BCReader.update_midmonth_data!(date, bcf_info)

        # nearest_idx = argmin(
        #     abs.(
        #         parse(FT, BCReader.datetime_to_strdate(date)) .-
        #         parse.(FT, BCReader.datetime_to_strdate.(bcf_info.all_dates[:]))
        #     ),
        # )
        # @test bcf_info.segment_idx[1] == bcf_info.segment_idx0[1] == nearest_idx

        # #  case 4: everything else
        # reset_bcf_info(bcf_info)
        # bcf_info.segment_idx[1] = bcf_info.segment_idx0[1] + Int(1)
        # date = bcf_info.all_dates[bcf_info.segment_idx[1]] - Dates.Day(1)

        # @test_throws ErrorException BCReader.update_midmonth_data!(date, bcf_info)
#     end
"""


# end

# delete testing directory and files
rm(regrid_dir; recursive = true, force = true)

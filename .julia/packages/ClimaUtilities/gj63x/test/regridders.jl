using Test

import ClimaUtilities
import ClimaUtilities: Regridders
import NCDatasets
import ClimaCore
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends

const context = ClimaComms.context()
ClimaComms.init(context)

include("TestTools.jl")

@testset "default_regridder_type" begin
    # Case 1: no regridder available
    @test_throws ErrorException Regridders.default_regridder_type()

    # Case 2: only TempestRegridder available
    import ClimaCoreTempestRemap
    @test Regridders.default_regridder_type() == :TempestRegridder

    # Case 3: TempestRegridder and InterpolationsRegridder both available
    import Interpolations
    @test Regridders.default_regridder_type() == :InterpolationsRegridder

    # Case 4: only TempestRegridder available
    # This case is not currently tested because we don't have a way to remove
    #  previously-loaded extensions.
end

@testset "InterpolationsRegridder incorrect dimensions" begin
    lon, lat, z =
        collect(-180.0:1:180.0), collect(-90.0:1:90), collect(0.0:1.0:100.0)
    dimensions3D = (lon, lat, z)
    dimensions3D_reversed = (lon, reverse(lat), reverse(z))
    dimensions2D = (lon, lat)
    dimensions2D_reversed = (lon, reverse(lat))
    size3D = (361, 181, 101)
    size2D = (361, 181)
    data_lat2D = zeros(size2D)
    data_lat2D_reversed = zeros(size2D)
    data_lon2D = zeros(size2D)
    data_lat3D = zeros(size3D)
    data_lat3D_reversed = zeros(size3D)
    data_lon3D = zeros(size3D)
    data_z3D = zeros(size3D)
    data_z3D_reversed = zeros(size3D)
    for i in 1:length(lon)
        data_lat2D[i, :] .= lat
        data_lat2D_reversed[i, :] .= reverse(lat)
    end
    for i in 1:length(lat)
        data_lon2D[:, i] .= lon
    end
    for i in 1:length(lon)
        for j in 1:length(z)
            data_lat3D_reversed[i, :, :] .= reverse(lat)
            data_lat3D[i, :, :] .= lat
        end
    end
    for i in 1:length(lat)
        for j in 1:length(z)
            data_lon3D[:, i, :] .= lon
        end
    end
    for i in 1:length(lon)
        for j in 1:length(lat)
            data_z3D_reversed[i, j, :] .= reverse(z)
            data_z3D[i, j, :] .= z
        end
    end
    spaces = make_spherical_space(Float64; context)
    horzspace = spaces.horizontal
    hv_center_space = spaces.hybrid


    # 3D space
    extrapolation_bc = (
        Interpolations.Periodic(),
        Interpolations.Flat(),
        Interpolations.Flat(),
    )
    # create one regirdder with no transformations to dimensions needed
    # create another regridder that reverses the second dimension
    reg_horz = Regridders.InterpolationsRegridder(horzspace)
    reg_horz_reversed = Regridders.InterpolationsRegridder(
        horzspace;
        dim_increasing = (true, false),
    )
    # check that `reg_horz_reversed` reverses lon and not lat as expected
    regridded_lat = Regridders.regrid(reg_horz, data_lat2D, dimensions2D)
    regridded_lon = Regridders.regrid(reg_horz, data_lon2D, dimensions2D)
    regridded_lat_reversed = Regridders.regrid(
        reg_horz_reversed,
        data_lat2D_reversed,
        dimensions2D_reversed,
    )
    regridded_lon_reversed =
        Regridders.regrid(reg_horz_reversed, data_lon2D, dimensions2D_reversed)
    @test regridded_lat_reversed == regridded_lat
    @test regridded_lon_reversed == regridded_lon

    # Create one regridder with no transformations to dimensions needed
    # Create another regridder that reverses the second and third dimensions
    reg_hv =
        Regridders.InterpolationsRegridder(hv_center_space; extrapolation_bc)
    regridded_lat = Regridders.regrid(reg_hv, data_lat3D, dimensions3D)
    regridded_lon = Regridders.regrid(reg_hv, data_lon3D, dimensions3D)
    regridded_z = Regridders.regrid(reg_hv, data_z3D, dimensions3D)
    dim_increasing = (true, false, false)
    reg_hv_reversed = Regridders.InterpolationsRegridder(
        hv_center_space;
        extrapolation_bc,
        dim_increasing,
    )
    regridded_lat_reversed = Regridders.regrid(
        reg_hv_reversed,
        data_lat3D_reversed,
        dimensions3D_reversed,
    )
    regridded_lon_reversed =
        Regridders.regrid(reg_hv_reversed, data_lon3D, dimensions3D_reversed)
    regridded_z_reversed = Regridders.regrid(
        reg_hv_reversed,
        data_z3D_reversed,
        dimensions3D_reversed,
    )
    @test regridded_lat_reversed == regridded_lat
    @test regridded_lon_reversed == regridded_lon
    @test regridded_z_reversed == regridded_z

    @test_throws "Dimensions must be monotonically increasing to use regrid!. Sort the dimensions first, or use regrid." Regridders.regrid!(
        zeros(axes(regridded_z_reversed)),
        reg_hv_reversed,
        data_lat3D,
        dimensions3D_reversed,
    )
end

@testset "InterpolationsRegridder" begin

    lon, lat, z =
        collect(-180.0:1:180.0), collect(-90.0:1:90), collect(0.0:1.0:100.0)
    dimensions2D = (lon, lat)
    dimensions3D = (lon, lat, z)
    size2D = (361, 181)
    size3D = (361, 181, 101)
    data_lat2D = zeros(size2D)
    data_lon2D = zeros(size2D)
    for i in 1:length(lon)
        data_lat2D[i, :] .= lat
    end
    for i in 1:length(lat)
        data_lon2D[:, i] .= lon
    end

    data_lat3D = zeros(size3D)
    data_lon3D = zeros(size3D)
    data_z3D = zeros(size3D)
    for i in 1:length(lon)
        for j in 1:length(z)
            data_lat3D[i, :, :] .= lat
        end
    end
    for i in 1:length(lat)
        for j in 1:length(z)
            data_lon3D[:, i, :] .= lon
        end
    end
    for i in 1:length(lon)
        for j in 1:length(lat)
            data_z3D[i, j, :] .= z
        end
    end

    for FT in (Float32, Float64)
        spaces = make_spherical_space(FT; context)
        horzspace = spaces.horizontal
        hv_center_space = spaces.hybrid

        reg_horz = Regridders.InterpolationsRegridder(horzspace)

        regridded_lat = Regridders.regrid(reg_horz, data_lat2D, dimensions2D)
        regridded_lon = Regridders.regrid(reg_horz, data_lon2D, dimensions2D)

        # now repeat the regrid with in place method, and check that the result is the same
        in_place_regridded_lat = zeros(axes(regridded_lat))
        in_place_regridded_lon = zeros(axes(regridded_lon))
        Regridders.regrid!(
            in_place_regridded_lat,
            reg_horz,
            FT.(data_lat2D),
            map(x -> FT.(x), dimensions2D),
        )
        Regridders.regrid!(
            in_place_regridded_lon,
            reg_horz,
            FT.(data_lon2D),
            map(x -> FT.(x), dimensions2D),
        )

        @test regridded_lat == in_place_regridded_lat
        @test regridded_lon == in_place_regridded_lon

        coordinates = ClimaCore.Fields.coordinate_field(horzspace)

        # Since lon uses periodic BCs but is not periodic itself, error is large at 180 degrees.
        # We check that the error is small at other indices, and that the values at -180 and 180 agree.
        function check_lon_error(coords_lon, lon_regrid)
            inds_normal = findall(!=(180), Array(parent(coords_lon)))
            err_lon = abs.(
                Array(parent(coords_lon))[inds_normal] .-
                Array(parent(lon_regrid))[inds_normal],
            )
            @test maximum(err_lon) < 1e-4

            inds_lon_180 = findall(==(180), Array(parent(coords_lon)))
            inds_lon_neg180 = findall(==(-180), Array(parent(coords_lon)))
            err_lon_180 = abs.(
                Array(parent(lon_regrid))[inds_lon_180] .-
                Array(parent(coords_lon))[inds_lon_neg180],
            )
            @test maximum(err_lon_180) < 1e-5
        end

        # Compute max err
        err_lat = abs.(coordinates.lat .- regridded_lat)
        @test maximum(err_lat) < 1e-5
        check_lon_error(coordinates.long, regridded_lon)

        # 3D space
        extrapolation_bc = (
            Interpolations.Periodic(),
            Interpolations.Flat(),
            Interpolations.Flat(),
        )

        reg_hv = Regridders.InterpolationsRegridder(
            hv_center_space;
            extrapolation_bc,
        )

        regridded_lat = Regridders.regrid(reg_hv, data_lat3D, dimensions3D)
        regridded_lon = Regridders.regrid(reg_hv, data_lon3D, dimensions3D)
        regridded_z = Regridders.regrid(reg_hv, data_z3D, dimensions3D)

        in_place_regridded_lat = zeros(axes(regridded_lat))
        in_place_regridded_lon = zeros(axes(regridded_lon))
        in_place_regridded_z = zeros(axes(regridded_z))
        Regridders.regrid!(
            in_place_regridded_lat,
            reg_hv,
            FT.(data_lat3D),
            map(x -> FT.(x), dimensions3D),
        )
        Regridders.regrid!(
            in_place_regridded_lon,
            reg_hv,
            FT.(data_lon3D),
            map(x -> FT.(x), dimensions3D),
        )
        Regridders.regrid!(
            in_place_regridded_z,
            reg_hv,
            FT.(data_z3D),
            map(x -> FT.(x), dimensions3D),
        )

        @test regridded_lat == in_place_regridded_lat
        @test regridded_lon == in_place_regridded_lon
        @test regridded_z == in_place_regridded_z

        coordinates = ClimaCore.Fields.coordinate_field(hv_center_space)

        # Compute max err
        err_lat = abs.(coordinates.lat .- regridded_lat)
        err_z = abs.(coordinates.z .- regridded_z)

        @test maximum(err_lat) < 1e-5
        @test maximum(err_z) < 1e-5
        check_lon_error(coordinates.long, regridded_lon)

        # 2D space with LatLongZ coordinates
        surface_space = ClimaCore.Spaces.level(hv_center_space, 1) # SpectralElementSpace2D with LatLongZPoint coordinates
        coordinates = ClimaCore.Fields.coordinate_field(surface_space)
        @assert !(
            surface_space isa ClimaCore.Spaces.ExtrudedFiniteDifferenceSpace
        )
        @assert eltype(coordinates) <: ClimaCore.Geometry.LatLongZPoint

        # 3D data
        reg_2d = Regridders.InterpolationsRegridder(surface_space)

        regridded_lat_3d = Regridders.regrid(reg_2d, data_lat3D, dimensions3D)
        regridded_lon_3d = Regridders.regrid(reg_2d, data_lon3D, dimensions3D)

        # Compute max err
        err_lat = abs.(coordinates.lat .- regridded_lat_3d)
        @test maximum(err_lat) < 1e-5
        check_lon_error(coordinates.long, regridded_lon_3d)

        # 2D data
        regridded_lat_2d = Regridders.regrid(reg_2d, data_lat2D, dimensions2D)
        regridded_lon_2d = Regridders.regrid(reg_2d, data_lon2D, dimensions2D)

        err_lat = abs.(coordinates.lat .- regridded_lat_2d)
        @test maximum(err_lat) < 1e-5
        check_lon_error(coordinates.long, regridded_lon_2d)
    end
end

@testset "Interpolation method" begin
    # Test constant interpolation method
    lon, lat, z =
        collect(0.0:1:360), collect(-90.0:1:90), collect(0.0:1.0:100.0)
    dimensions2D = (lon, lat)
    dimensions3D = (lon, lat, z)
    size2D = (361, 181)
    size3D = (361, 181, 101)
    ones_2Ddata = ones(size2D)
    ones_3Ddata = ones(size3D)

    for FT in (Float32, Float64)
        spaces = make_spherical_space(FT; context)
        horzspace = spaces.horizontal
        hv_center_space = spaces.hybrid

        reg_horz = Regridders.InterpolationsRegridder(
            horzspace,
            interpolation_method = Interpolations.Constant(),
        )
        reg_hv = Regridders.InterpolationsRegridder(
            hv_center_space,
            interpolation_method = Interpolations.Constant(),
        )

        regridded_2Ddata =
            Regridders.regrid(reg_horz, ones_2Ddata, dimensions2D)
        regridded_3Ddata = Regridders.regrid(reg_hv, ones_3Ddata, dimensions3D)

        @test all(x -> x == one(x), parent(regridded_2Ddata))
        @test all(x -> x == one(x), parent(regridded_3Ddata))
    end
end

@testset "InterpolationsRegridderXYZPoint" begin
    helem = (10, 10)
    Nq = 4
    zelem = 10
    FT = Float64 # test for Float64

    x, y, z = collect(0.0:1:5), collect(0.0:1:6), collect(0.0:1:7)

    dimensions3D = (x, y, z)
    size3D = (6, 7, 8)

    data_x3D = zeros(size3D)
    data_y3D = zeros(size3D)
    data_z3D = zeros(size3D)

    for i in 1:length(x)
        for j in 1:length(y)
            for k in 1:length(z)
                data_x3D[i, j, k] = x[i]
                data_y3D[i, j, k] = y[j]
                data_z3D[i, j, k] = z[k]
            end
        end
    end
    # create the box space with helper function
    spaces = make_box_space(Float64; context)
    # create interpolation object with InterpolationsRegridder for XYZPoint object
    reg_box = Regridders.InterpolationsRegridder(spaces)

    regridded_x = Regridders.regrid(reg_box, data_x3D, dimensions3D)

    regridded_y = Regridders.regrid(reg_box, data_y3D, dimensions3D)

    regridded_z = Regridders.regrid(reg_box, data_z3D, dimensions3D)
    # repeat the regrid with in place method, and check that the result is the same
    regridded_x_inplace = zeros(axes(regridded_x))
    regridded_y_inplace = zeros(axes(regridded_y))
    regridded_z_inplace = zeros(axes(regridded_z))

    Regridders.regrid!(
        regridded_x_inplace,
        reg_box,
        FT.(data_x3D),
        map(x -> FT.(x), dimensions3D),
    )
    Regridders.regrid!(
        regridded_y_inplace,
        reg_box,
        FT.(data_y3D),
        map(x -> FT.(x), dimensions3D),
    )
    Regridders.regrid!(
        regridded_z_inplace,
        reg_box,
        FT.(data_z3D),
        map(x -> FT.(x), dimensions3D),
    )

    @test regridded_x == regridded_x_inplace
    @test regridded_y == regridded_y_inplace
    @test regridded_z == regridded_z_inplace

    err_x = reg_box.coordinates.x .- regridded_x
    err_y = reg_box.coordinates.y .- regridded_y
    err_z = reg_box.coordinates.z .- regridded_z

    @test maximum(err_x) < 1e-5
    @test maximum(err_y) < 1e-5
    @test maximum(err_z) < 1e-5
end

@testset "TempestRegridder" begin
    for FT in (Float32, Float64)
        data_path = joinpath(@__DIR__, "test_data", "era5_1979_1.0x1.0_lai.nc")
        ds = NCDatasets.NCDataset(data_path, "r")
        original_max = maximum(ds["lai_lv"][:, :, 1])
        original_min = minimum(ds["lai_lv"][:, :, 1])
        test_time = ds["time"][1]
        close(ds)
        test_space = make_spherical_space(FT; context).horizontal
        regrid_dir = nothing
        if ClimaComms.iamroot(context)
            regrid_dir = mktempdir()
        end
        regrid_dir = ClimaComms.bcast(context, regrid_dir)
        ClimaComms.barrier(context)
        regridder = Regridders.TempestRegridder(
            test_space,
            "lai_lv",
            data_path;
            regrid_dir,
        )
        ClimaComms.barrier(context)
        regridded_field = Regridders.regrid(regridder, test_time)
        regridded_max = maximum(regridded_field)
        regridded_min = minimum(regridded_field)
        @test original_max >= regridded_max
        @test original_min <= regridded_min
    end
end

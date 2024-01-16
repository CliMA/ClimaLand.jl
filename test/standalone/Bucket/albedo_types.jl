using Test

using ClimaCore
using ClimaCore: Geometry, Meshes, Domains, Topologies, Spaces, Fields
using ClimaComms
import CLIMAParameters as CP

using Dates
using NCDatasets
using JLD2
using ClimaLand.Regridder: read_from_hdf5
using ClimaLand.FileReader:
    next_date_in_file, to_datetime, nans_to_zero, get_data_at_date
using ClimaLand.Bucket:
    BucketModel,
    BucketModelParameters,
    BulkAlbedoFunction,
    BulkAlbedoStatic,
    BulkAlbedoTemporal,
    bareground_albedo_dataset_path,
    cesm2_albedo_dataset_path,
    set_initial_parameter_field!,
    next_albedo
using ClimaLand.Domains: coordinates, Column, SphericalShell
using ClimaLand:
    initialize,
    make_update_aux,
    make_set_initial_cache,
    PrescribedAtmosphere,
    PrescribedRadiativeFluxes,
    TimeVaryingInput

# Bucket model parameters
import ClimaLand
include(joinpath(pkgdir(ClimaLand), "parameters", "create_parameters.jl"))

function create_domain_2d(FT)
    rad = FT(100)
    h = FT(3.5)
    ne = (2, 10)
    np = 2
    return SphericalShell(;
        radius = rad,
        depth = h,
        nelements = ne,
        npolynomial = np,
    )
end

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
# Use two separate regrid dirs to avoid duplicate filenames since files have same varname
regrid_dir_static = joinpath(pkgdir(ClimaLand), "test", "static")
regrid_dir_temporal = joinpath(pkgdir(ClimaLand), "test", "temporal")
isdir(regrid_dir_static) ? nothing : mkpath(regrid_dir_static)
isdir(regrid_dir_temporal) ? nothing : mkpath(regrid_dir_temporal)

@testset "Test set_initial_parameter_field for BulkAlbedoFunction, FT = $FT" begin
    # set up for function call
    α_sfc = (coord_point) -> sin(coord_point.lat + coord_point.long)
    α_snow = FT(0.8)
    albedo = BulkAlbedoFunction(α_snow, α_sfc)

    domain = create_domain_2d(FT)
    space = domain.space.surface
    p = (; bucket = (; α_sfc = Fields.zeros(space)))
    surface_coords = Fields.coordinate_field(space)

    set_initial_parameter_field!(albedo, p, surface_coords)

    # compare calculated result to manually applied albedo function
    @test p.bucket.α_sfc == α_sfc.(surface_coords)
end

@testset "Test next_albedo for BulkAlbedoFunction, FT = $FT" begin
    # set up each argument for function call
    α_sfc = (coord_point) -> sin(coord_point.lat + coord_point.long)
    α_snow = FT(0.8)
    albedo = BulkAlbedoFunction(α_snow, α_sfc)

    σS_c = FT(0.2)
    parameters = (; σS_c = σS_c)

    domain = create_domain_2d(FT)
    space = domain.space.surface
    surface_coords = Fields.coordinate_field(space)
    p = (; bucket = (; α_sfc = α_sfc.(surface_coords)))

    σS = FT(0.1)
    Y = (; bucket = (; σS = σS))

    # extract fields for manual calculation
    p_α_sfc = p.bucket.α_sfc
    Y_σS = Y.bucket.σS
    param_σS_c = parameters.σS_c
    a_α_snow = albedo.α_snow

    next_alb_manual = @. (
        (1 - Y_σS / (Y_σS + param_σS_c)) * p_α_sfc +
        Y_σS / (Y_σS + param_σS_c) * a_α_snow
    )

    @test next_albedo(albedo, parameters, Y, p, FT(0)) == next_alb_manual
end

# Add tests which use TempestRemap here -
# TempestRemap is not built on Windows because of NetCDF support limitations
# `BulkAlbedoStatic` and `BulkAlbedoTemporal` use TR via a call to
#   `hdwrite_regridfile_rll_to_cgll`
if !Sys.iswindows()
    @testset "Test next_albedo for BulkAlbedoStatic, FT = $FT" begin
        # setup for test
        domain = create_domain_2d(FT)
        space = domain.space.surface
        surface_coords = Fields.coordinate_field(space)

        # set up each argument for function call
        α_snow = FT(0.8)
        albedo = BulkAlbedoStatic{FT}(regrid_dir_static, space, α_snow = α_snow)

        σS_c = FT(0.2)
        parameters = (; σS_c = σS_c)

        p = (; bucket = (; α_sfc = Fields.zeros(space)))
        set_initial_parameter_field!(albedo, p, surface_coords)

        σS = FT(0.1)
        Y = (; bucket = (; σS = σS))

        # extract fields for manual calculation
        p_α_sfc = p.bucket.α_sfc
        Y_σS = Y.bucket.σS
        param_σS_c = parameters.σS_c
        a_α_snow = albedo.α_snow

        next_alb_manual = @. (
            (1 - Y_σS / (Y_σS + param_σS_c)) * p_α_sfc +
            Y_σS / (Y_σS + param_σS_c) * a_α_snow
        )

        @test next_albedo(albedo, parameters, Y, p, FT(0)) == next_alb_manual
    end

    @testset "Test set_initial_parameter_field for BulkAlbedoStatic, FT = $FT" begin
        # set up for function call
        domain = create_domain_2d(FT)
        space = domain.space.surface
        p = (; bucket = (; α_sfc = Fields.zeros(space)))
        surface_coords = Fields.coordinate_field(space)

        albedo = BulkAlbedoStatic{FT}(regrid_dir_static, space)
        set_initial_parameter_field!(albedo, p, surface_coords)

        # set up for manual data reading
        infile_path = bareground_albedo_dataset_path()
        varname = "sw_alb"

        data_manual = get_data_at_date(albedo.α_sfc, space, varname)

        @test p.bucket.α_sfc == data_manual
    end

    @testset "Test set_initial_parameter_field for BulkAlbedoTemporal, FT = $FT" begin
        # set up for function call
        t_start = Float64(0)
        domain = create_domain_2d(FT)
        space = domain.space.surface

        infile_path = cesm2_albedo_dataset_path()
        date_ref = to_datetime(NCDataset(infile_path, "r") do ds
            ds["time"][1]
        end)

        albedo = BulkAlbedoTemporal{FT}(
            regrid_dir_temporal,
            date_ref,
            t_start,
            space,
        )
        p = (; bucket = (; α_sfc = Fields.zeros(space)))
        surface_coords = Fields.coordinate_field(space)

        set_initial_parameter_field!(albedo, p, surface_coords)

        # set up for manual data reading
        infile_path = bareground_albedo_dataset_path()
        varname = "sw_alb"

        data_manual =
            get_data_at_date(albedo.albedo_info, space, varname, date_ref)

        @test nans_to_zero.(p.bucket.α_sfc) == nans_to_zero.(data_manual)
    end

    @testset "Test next_albedo for BulkAlbedoTemporal, FT = $FT" begin
        # set up each argument for function call
        domain = create_domain_2d(FT)
        space = domain.space.surface
        surface_coords = Fields.coordinate_field(space)

        infile_path = cesm2_albedo_dataset_path()
        date_ref = to_datetime(NCDataset(infile_path, "r") do ds
            ds["time"][1]
        end)
        t_start = Float64(0)

        albedo = BulkAlbedoTemporal{FT}(
            regrid_dir_temporal,
            date_ref,
            t_start,
            space,
        )

        Y = (; bucket = (; W = Fields.zeros(space)))
        p = (; bucket = (; α_sfc = Fields.zeros(space)))

        # initialize data fields
        set_initial_parameter_field!(albedo, p, surface_coords)

        # set up for manual data reading
        varname = "sw_alb"
        outfile_root = string(varname, "_cgll")
        file_dates = albedo.albedo_info.file_info.all_dates

        new_date = date_ref + Second(t_start)
        t_curr = t_start
        for i in 1:5
            @assert new_date == file_dates[i]

            # manually read in data at this date (not testing interpolation)
            field =
                get_data_at_date(albedo.albedo_info, space, varname, new_date)
            albedo_next = next_albedo(albedo, (;), Y, p, t_curr)
            replace_nan_missing!(albedo_next)
            replace_nan_missing!(field)
            @test albedo_next == field

            # Update manual date to match next date in file
            dt = Second(file_dates[i + 1] - file_dates[i])
            new_date += dt
            t_curr += dt.value
        end
    end

    @testset "Test BulkAlbedoStatic - albedo from map, FT = $FT" begin
        earth_param_set = create_lsm_parameters(FT)
        varname = "sw_alb"
        path = bareground_albedo_dataset_path()
        regrid_dirpath = joinpath(pkgdir(ClimaLand), "test/albedo_tmpfiles/")
        mkpath(regrid_dirpath)

        σS_c = FT(0.2)
        W_f = FT(0.15)
        z_0m = FT(1e-2)
        z_0b = FT(1e-3)
        κ_soil = FT(1.5)
        ρc_soil = FT(2e6)
        init_temp(z::FT, value::FT) where {FT} = FT(value)

        bucket_domains = [
            Column(; zlim = FT.((-100.0, 0.0)), nelements = 10),
            SphericalShell(;
                radius = FT(100.0),
                depth = FT(3.5),
                nelements = (2, 10),
                npolynomial = 2,
            ),
        ]

        for bucket_domain in bucket_domains
            if bucket_domain isa SphericalShell
                surface_space = bucket_domain.space.surface
                albedo_model =
                    BulkAlbedoStatic{FT}(regrid_dirpath, surface_space)

                # Radiation
                ref_time = DateTime(2005)
                SW_d = (t) -> 0.0
                LW_d = (t) -> 5.67e-8 * 280.0^4.0
                bucket_rad = PrescribedRadiativeFluxes(
                    FT,
                    TimeVaryingInput(SW_d),
                    TimeVaryingInput(LW_d),
                    ref_time,
                )
                # Atmos
                precip = (t) -> 0 # no precipitation
                T_atmos = (t) -> 280.0
                u_atmos = (t) -> 1.0
                q_atmos = (t) -> 0.0 # no atmos water
                h_atmos = FT(1e-8)
                P_atmos = (t) -> 101325
                bucket_atmos = PrescribedAtmosphere(
                    TimeVaryingInput(precip),
                    TimeVaryingInput(precip),
                    TimeVaryingInput(T_atmos),
                    TimeVaryingInput(u_atmos),
                    TimeVaryingInput(q_atmos),
                    TimeVaryingInput(P_atmos),
                    ref_time,
                    h_atmos,
                )
                τc = FT(1.0)
                bucket_parameters = BucketModelParameters(
                    κ_soil,
                    ρc_soil,
                    albedo_model,
                    σS_c,
                    W_f,
                    z_0m,
                    z_0b,
                    τc,
                    earth_param_set,
                )

                model = BucketModel(
                    parameters = bucket_parameters,
                    domain = bucket_domain,
                    atmosphere = bucket_atmos,
                    radiation = bucket_rad,
                )
                # Initial conditions with no moisture
                Y, p, coords = initialize(model)
                Y.bucket.T .= init_temp.(coords.subsurface.z, FT(280.0))
                Y.bucket.W .= 0.0
                Y.bucket.Ws .= 0.0
                Y.bucket.σS .= 0.0
                set_initial_cache! = make_set_initial_cache(model)
                set_initial_cache!(p, Y, FT(0.0))
                field = get_data_at_date(
                    albedo_model.α_sfc,
                    model.domain.space.surface,
                    varname,
                )
                @test p.bucket.α_sfc == field
            else
                surface_space = bucket_domain.space.surface
                @test_throws "Using an albedo map requires a global run." BulkAlbedoStatic{
                    FT,
                }(
                    regrid_dirpath,
                    surface_space,
                )
            end
        end
        rm(regrid_dirpath, recursive = true)
    end

    @testset "Test BulkAlbedoTemporal error with static map, FT = $FT" begin
        get_infile = bareground_albedo_dataset_path
        date_ref = Dates.DateTime(1900, 1, 1)
        t_start = Float64(0)

        # Test error for non-time varying data
        domain = create_domain_2d(FT)
        space = domain.space.surface

        @test_throws "Using a temporal albedo map requires data with time dimension." BulkAlbedoTemporal{
            FT,
        }(
            regrid_dir_temporal,
            date_ref,
            t_start,
            space,
            get_infile = get_infile,
        )
    end

    @testset "Test BulkAlbedoTemporal - albedo from map over time, FT = $FT" begin
        earth_param_set = create_lsm_parameters(FT)
        varname = "sw_alb"
        infile_path = cesm2_albedo_dataset_path()
        regrid_dirpath = joinpath(pkgdir(ClimaLand), "test/albedo_tmpfiles/")
        mkpath(regrid_dirpath)

        σS_c = FT(0.2)
        W_f = FT(0.15)
        z_0m = FT(1e-2)
        z_0b = FT(1e-3)
        κ_soil = FT(1.5)
        ρc_soil = FT(2e6)
        init_temp(z::FT, value::FT) where {FT} = FT(value)

        t_start = Float64(0)
        date_ref = to_datetime(NCDataset(infile_path, "r") do ds
            ds["time"][1]
        end)

        bucket_domains = [
            Column(; zlim = FT.((-100.0, 0.0)), nelements = 10),
            SphericalShell(;
                radius = FT(100.0),
                depth = FT(3.5),
                nelements = (2, 10),
                npolynomial = 2,
            ),
        ]

        for bucket_domain in bucket_domains
            space = bucket_domain.space.surface
            if bucket_domain isa SphericalShell
                albedo_model = BulkAlbedoTemporal{FT}(
                    regrid_dirpath,
                    date_ref,
                    t_start,
                    space,
                )
                # Radiation
                ref_time = DateTime(2005)
                SW_d = (t) -> 0
                LW_d = (t) -> 5.67e-8 * 280.0^4.0
                bucket_rad = PrescribedRadiativeFluxes(
                    FT,
                    TimeVaryingInput(SW_d),
                    TimeVaryingInput(LW_d),
                    ref_time,
                )
                # Atmos
                precip = (t) -> 0 # no precipitation
                T_atmos = (t) -> 280.0
                u_atmos = (t) -> 1.0
                q_atmos = (t) -> 0.0 # no atmos water
                h_atmos = FT(1e-8)
                P_atmos = (t) -> 101325
                ref_time = DateTime(2005)
                bucket_atmos = PrescribedAtmosphere(
                    TimeVaryingInput(precip),
                    TimeVaryingInput(precip),
                    TimeVaryingInput(T_atmos),
                    TimeVaryingInput(u_atmos),
                    TimeVaryingInput(q_atmos),
                    TimeVaryingInput(P_atmos),
                    ref_time,
                    h_atmos,
                )
                τc = FT(1.0)
                bucket_parameters = BucketModelParameters(
                    κ_soil,
                    ρc_soil,
                    albedo_model,
                    σS_c,
                    W_f,
                    z_0m,
                    z_0b,
                    τc,
                    earth_param_set,
                )

                model = BucketModel(
                    parameters = bucket_parameters,
                    domain = bucket_domain,
                    atmosphere = bucket_atmos,
                    radiation = bucket_rad,
                )
                # Initial conditions with no moisture
                Y, p, coords = initialize(model)
                Y.bucket.T .= init_temp.(coords.subsurface.z, FT(280.0))
                Y.bucket.W .= 0.0
                Y.bucket.Ws .= 0.0
                Y.bucket.σS .= 0.0
                set_initial_cache! = make_set_initial_cache(model)
                set_initial_cache!(p, Y, FT(0.0))
                data_manual = get_data_at_date(
                    albedo_model.albedo_info,
                    model.domain.space.surface,
                    varname,
                    date_ref,
                )
                # If there are any NaNs in the input data, replace them so we can compare results
                @test nans_to_zero.(p.bucket.α_sfc) ==
                      nans_to_zero.(data_manual)

                outfile_root = "temporal_data_cgll"
                file_dates = load(
                    joinpath(
                        regrid_dirpath,
                        outfile_root * "_" * varname * "_times.jld2",
                    ),
                    "times",
                )

                update_aux! = make_update_aux(model)
                new_date = date_ref + Second(t_start)
                t_curr = t_start
                for i in 1:5
                    @assert new_date == file_dates[i]

                    update_aux!(p, Y, t_curr)
                    data_manual = get_data_at_date(
                        albedo_model.albedo_info,
                        model.domain.space.surface,
                        varname,
                        file_dates[i],
                    )
                    @test nans_to_zero.(p.bucket.α_sfc) ≈
                          nans_to_zero.(data_manual)

                    # Update manual date to match next date in file
                    dt = Second(file_dates[i + 1] - file_dates[i])
                    new_date += dt
                    t_curr += dt.value
                end
            else
                @test_throws "Using an albedo map requires a global run." BulkAlbedoTemporal{
                    FT,
                }(
                    regrid_dirpath,
                    date_ref,
                    t_start,
                    space,
                )
            end
        end
        rm(regrid_dirpath, recursive = true)
    end
end

# Delete testing directory and files
rm(regrid_dir_static; recursive = true, force = true)
rm(regrid_dir_temporal; recursive = true, force = true)

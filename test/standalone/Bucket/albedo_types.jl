using Test

using ClimaCore
using ClimaCore: Geometry, Meshes, Domains, Topologies, Spaces, Fields
using ClimaComms
import ClimaParams as CP

using Dates
using NCDatasets
using ClimaLand.Regridder: read_from_hdf5
using ClimaLand.FileReader: next_date_in_file, to_datetime, get_data_at_date
using ClimaLand.Bucket:
    BucketModel,
    BucketModelParameters,
    PrescribedBaregroundAlbedo,
    PrescribedSurfaceAlbedo,
    bareground_albedo_dataset_path,
    cesm2_albedo_dataset_path,
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
import ClimaLand.Parameters as LP

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

@testset "Test next_albedo for PrescribedBaregroundAlbedo, from a function, FT = $FT" begin
    # set up each argument for function call
    α_bareground_func = (coord_point) -> sin(coord_point.lat + coord_point.long)
    α_snow = FT(0.8)
    domain = create_domain_2d(FT)
    space = domain.space.surface
    albedo = PrescribedBaregroundAlbedo{FT}(α_snow, α_bareground_func, space)

    σS_c = FT(0.2)
    parameters = (; σS_c = σS_c)

    surface_coords = Fields.coordinate_field(space)
    σS = FT(0.1)
    Y = (; bucket = (; σS = σS))
    p = (; bucket = (; α_sfc = Fields.zeros(space)))

    # extract fields for manual calculation
    α_bareground = α_bareground_func.(surface_coords)
    Y_σS = Y.bucket.σS
    param_σS_c = parameters.σS_c
    a_α_snow = albedo.α_snow

    next_alb_manual = @. (
        (1 - Y_σS / (Y_σS + param_σS_c)) * α_bareground +
        Y_σS / (Y_σS + param_σS_c) * a_α_snow
    )

    @test next_albedo(albedo, parameters, Y, p, FT(0)) == next_alb_manual
end

# Add tests which use TempestRemap here -
# TempestRemap is not built on Windows because of NetCDF support limitations
# Our prescribed albedo models  use TR via a call to
#   `hdwrite_regridfile_rll_to_cgll`
if !Sys.iswindows()
    @testset "Test next_albedo for PrescribedBaregroundAlbedo, from a map, FT = $FT" begin
        # setup for test
        domain = create_domain_2d(FT)
        space = domain.space.surface
        surface_coords = Fields.coordinate_field(space)

        # set up each argument for function call
        α_snow = FT(0.8)
        albedo =
            PrescribedBaregroundAlbedo{FT}(α_snow, regrid_dir_static, space)

        σS_c = FT(0.2)
        parameters = (; σS_c = σS_c)

        σS = FT(0.1)
        Y = (; bucket = (; σS = σS))
        p = (; bucket = (; α_sfc = Fields.zeros(space)))

        # extract fields for manual calculation
        Y_σS = Y.bucket.σS
        param_σS_c = parameters.σS_c
        a_α_snow = albedo.α_snow
        a_α_bareground = albedo.α_bareground
        @test eltype(a_α_bareground) == FT

        next_alb_manual = @. (
            (1 - Y_σS / (Y_σS + param_σS_c)) * a_α_bareground +
            Y_σS / (Y_σS + param_σS_c) * a_α_snow
        )

        @test next_albedo(albedo, parameters, Y, p, FT(0)) == next_alb_manual
    end

    @testset "Test next_albedo for PrescribedSurfaceAlbedo, FT = $FT" begin
        # set up each argument for function call
        domain = create_domain_2d(FT)
        space = domain.space.surface
        surface_coords = Fields.coordinate_field(space)

        infile_path = cesm2_albedo_dataset_path()
        date_ref = to_datetime(NCDataset(infile_path, "r") do ds
            ds["time"][1]
        end)
        t_start = Float64(0)

        albedo = PrescribedSurfaceAlbedo{FT}(
            regrid_dir_temporal,
            date_ref,
            t_start,
            space,
        )

        Y = (; bucket = (; W = Fields.zeros(space)))
        p = (; bucket = (; α_sfc = Fields.zeros(space)))

        # set up for manual data reading
        varname = "sw_alb"
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

    @testset "Test PrescribedBaregroundAlbedo - albedo from map, FT = $FT" begin
        earth_param_set = LP.LandParameters(FT)
        varname = "sw_alb"
        path = bareground_albedo_dataset_path()
        regrid_dirpath = joinpath(pkgdir(ClimaLand), "test/albedo_tmpfiles/")
        mkpath(regrid_dirpath)
        α_snow = FT(0.8)
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
                albedo = PrescribedBaregroundAlbedo{FT}(
                    α_snow,
                    regrid_dirpath,
                    surface_space,
                )

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
                bucket_parameters =
                    BucketModelParameters(FT; albedo, z_0m, z_0b, τc)

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

                # Read in data manually
                outfile_root = "static_data_cgll"
                comms_ctx = ClimaComms.context(surface_space)
                field = read_from_hdf5(
                    regrid_dirpath,
                    outfile_root,
                    Dates.DateTime(0), # dummy date
                    varname,
                    surface_space,
                )
                replace_nan_missing!(field)
                data_manual = field

                @test p.bucket.α_sfc == data_manual
            else
                surface_space = bucket_domain.space.surface
                @test_throws "Using an albedo map requires a global run." PrescribedBaregroundAlbedo{
                    FT,
                }(
                    α_snow,
                    regrid_dirpath,
                    surface_space,
                )
            end
        end
        rm(regrid_dirpath, recursive = true)
    end

    @testset "Test PrescribedSurfaceAlbedo error with static map, FT = $FT" begin
        get_infile = bareground_albedo_dataset_path
        date_ref = Dates.DateTime(1900, 1, 1)
        t_start = Float64(0)

        # Test error for non-time varying data
        domain = create_domain_2d(FT)
        space = domain.space.surface

        @test_throws "Using a temporal albedo map requires data with time dimension." PrescribedSurfaceAlbedo{
            FT,
        }(
            regrid_dir_temporal,
            date_ref,
            t_start,
            space,
            get_infile = get_infile,
        )
    end

    @testset "Test PrescribedSurfaceAlbedo - albedo from map over time, FT = $FT" begin
        earth_param_set = LP.LandParameters(FT)
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
        file_dates = to_datetime(NCDataset(infile_path, "r") do ds
            ds["time"][:]
        end)
        date_ref = file_dates[1]

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
                albedo_model = PrescribedSurfaceAlbedo{FT}(
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
                    FT;
                    albedo = albedo_model,
                    z_0m,
                    z_0b,
                    τc,
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
                replace_nan_missing!(p.bucket.α_sfc)
                replace_nan_missing!(data_manual)
                @test p.bucket.α_sfc == data_manual

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
                    replace_nan_missing!(p.bucket.α_sfc)
                    replace_nan_missing!(data_manual)
                    @test p.bucket.α_sfc == data_manual

                    # Update manual date to match next date in file
                    dt = Second(file_dates[i + 1] - file_dates[i])
                    new_date += dt
                    t_curr += dt.value
                end
            else
                @test_throws "Using an albedo map requires a global run." PrescribedSurfaceAlbedo{
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

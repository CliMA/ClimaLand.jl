using Test
import ClimaComms
ClimaComms.@import_required_backends

using ClimaUtilities.TimeManager
using ClimaUtilities.DataHandling
using ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
using ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
using ClimaCore
using ClimaCore: Geometry, Meshes, Domains, Topologies, Spaces, Fields
using ClimaComms
import NCDatasets, Interpolations
import ClimaParams as CP

using Dates
using NCDatasets
using ClimaLand.Bucket:
    BucketModel,
    BucketModelParameters,
    PrescribedBaregroundAlbedo,
    PrescribedSurfaceAlbedo,
    next_albedo!
using ClimaLand.Domains: coordinates, Column, SphericalShell
using ClimaLand:
    initialize,
    make_update_aux,
    make_set_initial_cache,
    PrescribedAtmosphere,
    PrescribedRadiativeFluxes


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

name_ds_var_list = [
    ("CESM2", ClimaLand.Artifacts.cesm2_albedo_dataset_path(), "sw_alb"),
    ("CERES", ClimaLand.Artifacts.ceres_albedo_dataset_path(), "sw_alb"),
    (
        "CERES CLR",
        ClimaLand.Artifacts.ceres_albedo_dataset_path(),
        "sw_alb_clr",
    ),
]

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

    next_albedo!(p.bucket.α_sfc, albedo, parameters, Y, p, FT(0))

    @test p.bucket.α_sfc == next_alb_manual
end

@testset "Test next_albedo for PrescribedBaregroundAlbedo, from a map, FT = $FT" begin
    # setup for test
    domain = create_domain_2d(FT)
    space = domain.space.surface
    surface_coords = Fields.coordinate_field(space)

    # set up each argument for function call
    α_snow = FT(0.8)
    albedo = PrescribedBaregroundAlbedo{FT}(α_snow, space)

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

    next_albedo!(p.bucket.α_sfc, albedo, parameters, Y, p, FT(0))
    @test p.bucket.α_sfc == next_alb_manual
end

for (name, infile_path, varname) in name_ds_var_list
    @testset "Test next_albedo for PrescribedSurfaceAlbedo, FT = $FT" begin
        # set up each argument for function call
        domain = create_domain_2d(FT)
        space = domain.space.surface
        surface_coords = Fields.coordinate_field(space)

        start_date_noleap = NCDataset(infile_path, "r") do ds
            ds["time"][1]
        end
        # Converting from NoLeap
        start_date = Dates.DateTime(
            Dates.year(start_date_noleap),
            Dates.month(start_date_noleap),
            Dates.day(start_date_noleap),
            Dates.hour(start_date_noleap),
            Dates.minute(start_date_noleap),
            Dates.second(start_date_noleap),
            Dates.millisecond(start_date_noleap),
        )
        t_start = Float64(0)

        albedo = PrescribedSurfaceAlbedo{FT}(
            start_date,
            space;
            albedo_file_path = infile_path,
            varname = varname,
        )

        Y = (; bucket = (; W = Fields.zeros(space)))
        p = (; bucket = (; α_sfc = Fields.zeros(space)))

        # set up for manual data reading
        file_dates = DataHandling.available_dates(albedo.albedo.data_handler)

        new_date = start_date + Second(t_start)
        t_curr = t_start
        for i in 1:5
            @assert new_date == file_dates[i]

            # manually read in data at this date (not testing interpolation)
            field = DataHandling.regridded_snapshot(
                albedo.albedo.data_handler,
                new_date,
            )

            next_albedo!(p.bucket.α_sfc, albedo, (;), Y, p, t_curr)
            @test p.bucket.α_sfc == field

            # Update manual date to match next date in file
            dt = Second(file_dates[i + 1] - file_dates[i])
            new_date += dt
            t_curr += dt.value
        end
    end
end

@testset "Test PrescribedBaregroundAlbedo - albedo from map, FT = $FT" begin
    earth_param_set = LP.LandParameters(FT)
    varname = "sw_alb"
    path = ClimaLand.Artifacts.bareground_albedo_dataset_path()
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
            albedo = PrescribedBaregroundAlbedo{FT}(α_snow, surface_space)

            # Radiation
            start_date = DateTime(2005, 1, 15, 12)
            SW_d = (t) -> 0.0
            LW_d = (t) -> 5.67e-8 * 280.0^4.0
            bucket_rad = PrescribedRadiativeFluxes(
                FT,
                TimeVaryingInput(SW_d),
                TimeVaryingInput(LW_d),
                start_date,
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
                start_date,
                h_atmos,
                earth_param_set,
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
            @test ClimaComms.context(model) == ClimaComms.context()
            @test ClimaComms.device(model) == ClimaComms.device()

            # Initial conditions with no moisture
            Y, p, coords = initialize(model)
            Y.bucket.T .= init_temp.(coords.subsurface.z, FT(280.0))
            Y.bucket.W .= 0.0
            Y.bucket.Ws .= 0.0
            Y.bucket.σS .= 0.0
            set_initial_cache! = make_set_initial_cache(model)
            set_initial_cache!(p, Y, FT(0.0))

            # Read in data manually
            α_bareground = SpaceVaryingInput(path, varname, surface_space)
            @test p.bucket.α_sfc == α_bareground
        else
            surface_space = bucket_domain.space.surface
            @test_throws "Using an albedo map requires a global run." PrescribedBaregroundAlbedo{
                FT,
            }(
                α_snow,
                surface_space,
            )
        end
    end
end

for (name, infile_path, varname) in name_ds_var_list
    @testset "Test PrescribedSurfaceAlbedo - albedo from map over time, FT = $FT" begin
        earth_param_set = LP.LandParameters(FT)

        σS_c = FT(0.2)
        W_f = FT(0.15)
        z_0m = FT(1e-2)
        z_0b = FT(1e-3)
        κ_soil = FT(1.5)
        ρc_soil = FT(2e6)
        init_temp(z::FT, value::FT) where {FT} = FT(value)

        t_start = Float64(0)
        file_dates_noleap = NCDataset(infile_path, "r") do ds
            ds["time"][:]
        end
        to_datetime(date) = Dates.DateTime(
            Dates.year(date),
            Dates.month(date),
            Dates.day(date),
            Dates.hour(date),
            Dates.minute(date),
            Dates.second(date),
            Dates.millisecond(date),
        )
        file_dates = to_datetime.(file_dates_noleap)
        start_date_in_file = file_dates[1]

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
                    start_date_in_file,
                    space;
                    albedo_file_path = infile_path,
                    varname = varname,
                )
                # Radiation
                start_date = DateTime(2005, 1, 15, 12)
                SW_d = (t) -> 0
                LW_d = (t) -> 5.67e-8 * 280.0^4.0
                bucket_rad = PrescribedRadiativeFluxes(
                    FT,
                    TimeVaryingInput(SW_d),
                    TimeVaryingInput(LW_d),
                    start_date,
                )
                # Atmos
                precip = (t) -> 0 # no precipitation
                T_atmos = (t) -> 280.0
                u_atmos = (t) -> 1.0
                q_atmos = (t) -> 0.0 # no atmos water
                h_atmos = FT(1e-8)
                P_atmos = (t) -> 101325
                start_date = DateTime(2005, 1, 15, 12)
                bucket_atmos = PrescribedAtmosphere(
                    TimeVaryingInput(precip),
                    TimeVaryingInput(precip),
                    TimeVaryingInput(T_atmos),
                    TimeVaryingInput(u_atmos),
                    TimeVaryingInput(q_atmos),
                    TimeVaryingInput(P_atmos),
                    start_date,
                    h_atmos,
                    earth_param_set,
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
                data_manual = DataHandling.regridded_snapshot(
                    albedo_model.albedo.data_handler,
                    start_date_in_file,
                )

                @test p.bucket.α_sfc == data_manual

                update_aux! = make_update_aux(model)
                new_date = start_date_in_file + Second(t_start)
                t_curr = t_start
                for i in 1:5
                    @assert new_date == file_dates[i]

                    update_aux!(p, Y, t_curr)
                    data_manual = DataHandling.regridded_snapshot(
                        albedo_model.albedo.data_handler,
                        new_date,
                    )

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
                    start_date_in_file,
                    space,
                )
            end
        end
    end
end

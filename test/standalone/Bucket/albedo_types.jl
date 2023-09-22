using Test

using ClimaCore
using ClimaCore: Geometry, Meshes, Domains, Topologies, Spaces, Fields
using ClimaComms
import CLIMAParameters as CP
using Insolation
using Dates
using NCDatasets
using JLD2
using ClimaLSM.Regridder: regrid_netcdf_to_field, nans_to_zero, read_from_hdf5
using ClimaLSM.FileReader: next_date_in_file, to_datetime
using ClimaLSM.Bucket:
    BucketModel,
    BucketModelParameters,
    BulkAlbedoFunction,
    BulkAlbedoStatic,
    BulkAlbedoTemporal,
    bareground_albedo_dataset_path,
    cesm2_albedo_dataset_path,
    set_initial_parameter_field!,
    next_albedo
using ClimaLSM.Domains: coordinates, Column, SphericalShell
using ClimaLSM:
    initialize,
    make_update_aux,
    make_set_initial_aux_state,
    PrescribedAtmosphere,
    PrescribedRadiativeFluxes

# Bucket model parameters
import ClimaLSM
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

# Use two separate regrid dirs to avoid duplicate filenames since files have same varname
regrid_dir_static = joinpath(pkgdir(ClimaLSM), "test", "static")
regrid_dir_temporal = joinpath(pkgdir(ClimaLSM), "test", "temporal")
isdir(regrid_dir_static) ? nothing : mkpath(regrid_dir_static)
isdir(regrid_dir_temporal) ? nothing : mkpath(regrid_dir_temporal)

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

@testset "Test set_initial_parameter_field for BulkAlbedoFunction" begin
    FT = Float32
    # set up for function call
    α_sfc = (coord_point) -> sin(coord_point.lat + coord_point.long)
    α_snow = FT(0.8)
    albedo = BulkAlbedoFunction{FT}(α_snow, α_sfc)

    domain = create_domain_2d(FT)
    space = domain.space.surface
    p = (; bucket = (; α_sfc = Fields.zeros(space)))
    surface_coords = Fields.coordinate_field(space)

    set_initial_parameter_field!(albedo, p, surface_coords)

    # compare calculated result to manually applied albedo function
    @test p.bucket.α_sfc == α_sfc.(surface_coords)
end

@testset "Test set_initial_parameter_field for BulkAlbedoStatic" begin
    FT = Float32
    # set up for function call
    regrid_dir_static = joinpath(pkgdir(ClimaLSM), "test", "static")
    albedo = BulkAlbedoStatic{FT}(regrid_dir_static)
    domain = create_domain_2d(FT)
    space = domain.space.surface
    p = (; bucket = (; α_sfc = Fields.zeros(space)))
    surface_coords = Fields.coordinate_field(space)

    set_initial_parameter_field!(albedo, p, surface_coords)

    # set up for manual data reading
    comms = ClimaComms.SingletonCommsContext()
    input_file = bareground_albedo_dataset_path()
    varname = "sw_alb"

    data_manual = regrid_netcdf_to_field(
        FT,
        regrid_dir_static,
        comms,
        input_file,
        varname,
        space,
    )

    @test p.bucket.α_sfc == data_manual
end

@testset "Test set_initial_parameter_field for BulkAlbedoTemporal" begin
    FT = Float32
    # set up for function call
    regrid_dir_temporal = joinpath(pkgdir(ClimaLSM), "test", "temporal")
    t_start = FT(0)
    domain = create_domain_2d(FT)
    space = domain.space.surface

    input_file = cesm2_albedo_dataset_path()
    date_ref = to_datetime(NCDataset(input_file, "r") do ds
        ds["time"][1]
    end)

    albedo =
        BulkAlbedoTemporal{FT}(regrid_dir_temporal, date_ref, t_start, space)
    p = (; bucket = (; α_sfc = Fields.zeros(space)))
    surface_coords = Fields.coordinate_field(space)

    set_initial_parameter_field!(albedo, p, surface_coords)

    # set up for manual data reading
    comms = ClimaComms.SingletonCommsContext()
    input_file = bareground_albedo_dataset_path()
    varname = "sw_alb"

    data_manual = regrid_netcdf_to_field(
        FT,
        regrid_dir_temporal,
        comms,
        input_file,
        varname,
        space,
    )

    @test nans_to_zero.(p.bucket.α_sfc) == nans_to_zero.(data_manual)
end

@testset "Test next_albedo for BulkAlbedoFunction" begin
    FT = Float32
    # set up each argument for function call
    α_sfc = (coord_point) -> sin(coord_point.lat + coord_point.long)
    α_snow = FT(0.8)
    albedo = BulkAlbedoFunction{FT}(α_snow, α_sfc)

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

@testset "Test next_albedo for BulkAlbedoStatic" begin
    FT = Float32
    # set up each argument for function call
    α_snow = FT(0.8)
    albedo = BulkAlbedoStatic{FT}(regrid_dir_static, α_snow = α_snow)

    σS_c = FT(0.2)
    parameters = (; σS_c = σS_c)

    domain = create_domain_2d(FT)
    space = domain.space.surface
    surface_coords = Fields.coordinate_field(space)
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

@testset "Test next_albedo for BulkAlbedoTemporal" begin
    FT = Float32
    # set up each argument for function call
    domain = create_domain_2d(FT)
    space = domain.space.surface
    surface_coords = Fields.coordinate_field(space)

    input_file = cesm2_albedo_dataset_path()
    date_ref = to_datetime(NCDataset(input_file, "r") do ds
        ds["time"][1]
    end)
    t_start = FT(0)

    albedo =
        BulkAlbedoTemporal{FT}(regrid_dir_temporal, date_ref, t_start, space)

    Y = (; bucket = (; W = Fields.zeros(space)))
    p = (; bucket = (; α_sfc = Fields.zeros(space)))

    # initialize data fields
    set_initial_parameter_field!(albedo, p, surface_coords)

    # set up for manual data reading
    varname = "sw_alb"
    outfile_root = string(varname, "_cgll")
    file_dates = load(
        joinpath(regrid_dir_temporal, outfile_root * "_times.jld2"),
        "times",
    )
    comms = ClimaComms.SingletonCommsContext()

    new_date = date_ref + Second(t_start)
    t_curr = t_start
    for i in 1:5
        @assert new_date == file_dates[i]

        # manually read in data at this date (not testing interpolation)
        field = regrid_netcdf_to_field(
            FT,
            regrid_dir_temporal,
            comms,
            input_file,
            varname,
            space,
            date_idx = i,
        )
        @test nans_to_zero.(next_albedo(albedo, (;), Y, (;), t_curr)) ≈ field

        # Update manual date to match next date in file
        dt = Second(file_dates[i + 1] - file_dates[i])
        new_date += dt
        t_curr += dt.value
    end
end

@testset "Test BulkAlbedoStatic - albedo from map" begin
    FT = Float32
    earth_param_set = create_lsm_parameters(FT)
    varname = "sw_alb"
    path = bareground_albedo_dataset_path()
    comms = ClimaComms.SingletonCommsContext()
    regrid_dirpath = joinpath(pkgdir(ClimaLSM), "test/albedo_tmpfiles/")
    mkpath(regrid_dirpath)
    albedo_model = BulkAlbedoStatic{FT}(regrid_dirpath)

    σS_c = FT(0.2)
    W_f = FT(0.15)
    z_0m = FT(1e-2)
    z_0b = FT(1e-3)
    κ_soil = FT(1.5)
    ρc_soil = FT(2e6)
    init_temp(z::FT, value::FT) where {FT} = FT(value)

    bucket_domains = [
        Column(; zlim = (-100.0, 0.0), nelements = 10),
        SphericalShell(;
            radius = FT(100.0),
            depth = FT(3.5),
            nelements = (2, 10),
            npolynomial = 2,
        ),
    ]
    orbital_data = Insolation.OrbitalData()


    for bucket_domain in bucket_domains
        # Radiation
        ref_time = DateTime(2005)
        SW_d = (t) -> eltype(t)(0.0)
        LW_d = (t) -> eltype(t)(5.67e-8 * 280.0^4.0)
        bucket_rad =
            PrescribedRadiativeFluxes(FT, SW_d, LW_d, ref_time; orbital_data)
        # Atmos
        precip = (t) -> eltype(t)(0) # no precipitation
        T_atmos = (t) -> eltype(t)(280.0)
        u_atmos = (t) -> eltype(t)(1.0)
        q_atmos = (t) -> eltype(t)(0.0) # no atmos water
        h_atmos = FT(1e-8)
        P_atmos = (t) -> eltype(t)(101325)
        bucket_atmos = PrescribedAtmosphere(
            precip,
            precip,
            T_atmos,
            u_atmos,
            q_atmos,
            P_atmos,
            ref_time,
            h_atmos,
        )
        Δt = FT(1.0)
        bucket_parameters = BucketModelParameters(
            κ_soil,
            ρc_soil,
            albedo_model,
            σS_c,
            W_f,
            z_0m,
            z_0b,
            Δt,
            earth_param_set,
        )
        if bucket_domain isa SphericalShell
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
            set_initial_aux_state! = make_set_initial_aux_state(model)
            set_initial_aux_state!(p, Y, FT(0.0))
            field = regrid_netcdf_to_field(
                FT,
                regrid_dirpath,
                comms,
                path,
                varname,
                model.domain.space.surface,
            )
            @test p.bucket.α_sfc == field
        else
            let err = nothing
                try
                    model = BucketModel(
                        parameters = bucket_parameters,
                        domain = bucket_domain,
                        atmosphere = bucket_atmos,
                        radiation = bucket_rad,
                    )
                catch err
                end

                @test err isa Exception
                @test sprint(showerror, err) ==
                      "Using an albedo map requires a global run."
            end
        end
    end
    rm(regrid_dirpath, recursive = true)
end

@testset "Test BulkAlbedoTemporal error with static map" begin
    FT = Float32
    regrid_dirpath = ""
    input_file = bareground_albedo_dataset_path()
    date_ref = Dates.DateTime(1900, 1, 1)
    t_start = FT(0)
    domain = create_domain_2d(FT)
    space = domain.space.surface

    let err = nothing
        try
            BulkAlbedoTemporal{FT}(
                regrid_dirpath,
                date_ref,
                t_start,
                space,
                input_file = input_file,
            )
        catch err
        end

        @test err isa Exception
        @test sprint(showerror, err) ==
              "Using a temporal albedo map requires data with time dimension."
    end

end

# Note: this test implicitly tests `FileReader.interpolate_data` behavior
@testset "Test BulkAlbedoTemporal - albedo from map over time" begin
    FT = Float32
    earth_param_set = create_lsm_parameters(FT)
    varname = "sw_alb"
    input_file = cesm2_albedo_dataset_path()
    comms = ClimaComms.SingletonCommsContext()
    regrid_dirpath = joinpath(pkgdir(ClimaLSM), "test/albedo_tmpfiles/")
    mkpath(regrid_dirpath)

    σS_c = FT(0.2)
    W_f = FT(0.15)
    z_0m = FT(1e-2)
    z_0b = FT(1e-3)
    κ_soil = FT(1.5)
    ρc_soil = FT(2e6)
    init_temp(z::FT, value::FT) where {FT} = FT(value)

    t_start = FT(0)
    date_ref = to_datetime(NCDataset(input_file, "r") do ds
        ds["time"][1]
    end)

    bucket_domains = [
        Column(; zlim = (-100.0, 0.0), nelements = 10),
        SphericalShell(;
            radius = FT(100.0),
            depth = FT(3.5),
            nelements = (2, 10),
            npolynomial = 2,
        ),
    ]
    orbital_data = Insolation.OrbitalData()

    for bucket_domain in bucket_domains
        space = bucket_domain.space.surface
        if bucket_domain isa SphericalShell
            albedo_model =
                BulkAlbedoTemporal{FT}(regrid_dirpath, date_ref, t_start, space)
            # Radiation
            ref_time = DateTime(2005)
            SW_d = (t) -> eltype(t)(0.0)
            LW_d = (t) -> eltype(t)(5.67e-8 * 280.0^4.0)
            bucket_rad = PrescribedRadiativeFluxes(
                FT,
                SW_d,
                LW_d,
                ref_time;
                orbital_data,
            )
            # Atmos
            precip = (t) -> eltype(t)(0) # no precipitation
            T_atmos = (t) -> eltype(t)(280.0)
            u_atmos = (t) -> eltype(t)(1.0)
            q_atmos = (t) -> eltype(t)(0.0) # no atmos water
            h_atmos = FT(1e-8)
            P_atmos = (t) -> eltype(t)(101325)
            ref_time = DateTime(2005)
            bucket_atmos = PrescribedAtmosphere(
                precip,
                precip,
                T_atmos,
                u_atmos,
                q_atmos,
                P_atmos,
                ref_time,
                h_atmos,
            )
            Δt = FT(1.0)
            bucket_parameters = BucketModelParameters(
                κ_soil,
                ρc_soil,
                albedo_model,
                σS_c,
                W_f,
                z_0m,
                z_0b,
                Δt,
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
            set_initial_aux_state! = make_set_initial_aux_state(model)
            set_initial_aux_state!(p, Y, FT(0.0))
            field = regrid_netcdf_to_field(
                FT,
                regrid_dirpath,
                comms,
                input_file,
                varname,
                model.domain.space.surface,
                date_idx = 1,
            )
            # If there are any NaNs in the input data, replace them so we can compare results
            @test nans_to_zero.(p.bucket.α_sfc) == field

            outfile_root = string(varname, "_cgll")
            file_dates = load(
                joinpath(regrid_dirpath, outfile_root * "_times.jld2"),
                "times",
            )

            update_aux! = make_update_aux(model)
            new_date = date_ref + Second(t_start)
            t_curr = t_start
            for i in 1:5
                @assert new_date == file_dates[i]

                update_aux!(p, Y, t_curr)
                field = regrid_netcdf_to_field(
                    FT,
                    regrid_dirpath,
                    comms,
                    input_file,
                    varname,
                    model.domain.space.surface,
                    date_idx = i,
                )
                @test nans_to_zero.(p.bucket.α_sfc) ≈ field

                # Update manual date to match next date in file
                dt = Second(file_dates[i + 1] - file_dates[i])
                new_date += dt
                t_curr += dt.value
            end
        else
            let err = nothing
                try
                    albedo_model = BulkAlbedoTemporal{FT}(
                        regrid_dirpath,
                        date_ref,
                        t_start,
                        space,
                    )
                catch err
                end

                @test err isa Exception
                @test sprint(showerror, err) ==
                      "Using an albedo map requires a global run."
            end
        end
    end
    rm(regrid_dirpath, recursive = true)
end

# Delete testing directory and files
rm(regrid_dir_static; recursive = true, force = true)
rm(regrid_dir_temporal; recursive = true, force = true)

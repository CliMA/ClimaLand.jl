using Test

using ClimaCore
using ClimaComms
import CLIMAParameters as CP
using Insolation
using Dates
using ClimaLSM.Regridder:
    MapInfo, regrid_netcdf_to_field, nans_to_zero, read_from_hdf5
using ClimaLSM.BCReader: next_date_in_file
using ClimaLSM.Bucket:
    BucketModel,
    BucketModelParameters,
    BulkAlbedoStatic,
    BulkAlbedoTemporal,
    bareground_albedo_dataset_path,
    cesm2_albedo_dataset_path
using ClimaLSM.Domains:
    coordinates, LSMSingleColumnDomain, LSMSphericalShellDomain
using ClimaLSM:
    initialize,
    make_update_aux,
    make_set_initial_aux_state,
    PrescribedAtmosphere,
    PrescribedRadiativeFluxes



# Bucket model parameters
import ClimaLSM
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

@testset "Test BulkAlbedoStatic - albedo from map" begin
    FT = Float32
    earth_param_set = create_lsm_parameters(FT)
    varname = "sw_alb"
    path = bareground_albedo_dataset_path()
    comms = ClimaComms.SingletonCommsContext()
    regrid_dirpath = joinpath(pkgdir(ClimaLSM), "test/Bucket/albedo_tmpfiles/")
    mkpath(regrid_dirpath)
    albedo_model = BulkAlbedoStatic{FT}(regrid_dirpath, path)

    σS_c = FT(0.2)
    W_f = FT(0.15)
    z_0m = FT(1e-2)
    z_0b = FT(1e-3)
    κ_soil = FT(1.5)
    ρc_soil = FT(2e6)
    init_temp(z::FT, value::FT) where {FT} = FT(value)

    bucket_domains = [
        LSMSingleColumnDomain(; zlim = (-100.0, 0.0), nelements = 10),
        LSMSphericalShellDomain(;
            radius = FT(100.0),
            height = FT(3.5),
            nelements = (2, 10),
            npolynomial = 2,
        ),
    ]
    orbital_data = Insolation.OrbitalData()


    for bucket_domain in bucket_domains
        # Radiation
        SW_d = (t) -> eltype(t)(0.0)
        LW_d = (t) -> eltype(t)(5.67e-8 * 280.0^4.0)
        bucket_rad = PrescribedRadiativeFluxes(FT, SW_d, LW_d; orbital_data)
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
        if bucket_domain isa LSMSphericalShellDomain
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
                model.domain.surface.space,
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
    path = bareground_albedo_dataset_path()
    surface = nothing
    err = nothing

    try BulkAlbedoTemporal{FT}(regrid_dirpath, path, surface)
    catch err
    end

    @test err isa Exception
    @test sprint(showerror, err) ==
            "Using a temporal albedo map requires data with time dimension."

end

@testset "Test BulkAlbedoTemporal - albedo from map over time" begin
    # function test_albedotemporal()
    FT = Float32
    earth_param_set = create_lsm_parameters(FT)
    varname = "sw_alb"
    path = cesm2_albedo_dataset_path()
    comms = ClimaComms.SingletonCommsContext()
    regrid_dirpath = joinpath(pkgdir(ClimaLSM), "test/Bucket/albedo_tmpfiles/")
    mkpath(regrid_dirpath)

    σS_c = FT(0.2)
    W_f = FT(0.15)
    z_0m = FT(1e-2)
    z_0b = FT(1e-3)
    κ_soil = FT(1.5)
    ρc_soil = FT(2e6)
    init_temp(z::FT, value::FT) where {FT} = FT(value)

    bucket_domains = [
        LSMSingleColumnDomain(; zlim = (-100.0, 0.0), nelements = 10),
        LSMSphericalShellDomain(;
            radius = FT(100.0),
            height = FT(3.5),
            nelements = (2, 10),
            npolynomial = 2,
        ),
    ]
    orbital_data = Insolation.OrbitalData()


    for bucket_domain in bucket_domains
        surface = bucket_domain.surface
        if bucket_domain isa LSMSphericalShellDomain
            albedo_model = BulkAlbedoTemporal{FT}(regrid_dirpath, path, surface)
            # Radiation
            SW_d = (t) -> eltype(t)(0.0)
            LW_d = (t) -> eltype(t)(5.67e-8 * 280.0^4.0)
            bucket_rad = PrescribedRadiativeFluxes(FT, SW_d, LW_d; orbital_data)
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
                path,
                varname,
                model.domain.surface.space,
            )
            # If there are any NaNs in the input data, replace them so we can compare results
            @test nans_to_zero.(p.bucket.α_sfc) == field
            # Update albedo over time
            t_start = FT(0)
            update_aux! = make_update_aux(model)
            update_aux!(p, Y, t_start)
            @test p.bucket.α_sfc != field

            # TODO debug this test - step to next date of file, check that α_sfc matches NCDataset file
            # next_date = next_date_in_file(albedo_model.albedo_info)
            # t_next = FT(Second(next_date - albedo_model.albedo_info.all_dates[1]).value)
            # update_aux!(p, Y, t_next)
            # next_field = read_from_hdf5(
            #                 regrid_dirpath,
            #                 varname * "cgll",
            #                 next_date,
            #                 varname,
            #                 comms,
            #             )
        else
            let err = nothing
                try
                    albedo_model =
                        BulkAlbedoTemporal{FT}(regrid_dirpath, path, surface)
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

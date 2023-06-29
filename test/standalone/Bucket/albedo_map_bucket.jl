using Test

using ClimaCore
using ClimaComms
import CLIMAParameters as CP
using Insolation
using ClimaLSM.Regridder: MapInfo, regrid_netcdf_to_field
using ClimaLSM.Bucket:
    BucketModel,
    BucketModelParameters,
    BulkAlbedoMap,
    bareground_albedo_dataset_path
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

@testset "Initialize aux state with map" begin
    FT = Float32
    earth_param_set = create_lsm_parameters(FT)
    varname = "sw_alb"
    path = bareground_albedo_dataset_path()
    comms = ClimaComms.SingletonCommsContext()
    regrid_dirpath =
        joinpath(pkgdir(ClimaLSM), "test/standalone/Bucket/albedo_tmpfiles")
    albedo_model = BulkAlbedoMap{FT}(regrid_dirpath)

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
        mkpath(regrid_dirpath)
        rm(regrid_dirpath, recursive = true)
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

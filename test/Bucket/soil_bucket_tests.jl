using Test

using Statistics
using ClimaCore
import CLIMAParameters as CP
using ClimaLSM.Bucket: BucketModel, BucketModelParameters, BulkAlbedoFunction
using ClimaLSM.Domains:
    coordinates,
    LSMSingleColumnDomain,
    LSMMultiColumnDomain,
    LSMSphericalShellDomain
using ClimaLSM:
    initialize,
    make_update_aux,
    make_ode_function,
    make_set_initial_aux_state,
    PrescribedAtmosphere,
    PrescribedRadiativeFluxes


FT = Float64

# Bucket model parameters
import ClimaLSM
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))
earth_param_set = create_lsm_parameters(FT)

α_sfc = (coordinate_point) -> FT(0.2) # surface albedo, spatially constant
α_snow = FT(0.8) # snow albedo
albedo = BulkAlbedoFunction{FT}(α_snow, α_sfc)
σS_c = FT(0.2)
W_f = FT(0.15)
z_0m = FT(1e-2)
z_0b = FT(1e-3)
κ_soil = FT(1.5)
ρc_soil = FT(2e6)

# Model domain
bucket_domains = [
    LSMSingleColumnDomain(; zlim = (-100.0, 0.0), nelements = 10),
    LSMMultiColumnDomain(;
        xlim = (-1.0, 0.0),
        ylim = (-1.0, 0.0),
        zlim = (-100.0, 0.0),
        nelements = (2, 2, 10),
        npolynomial = 1,
        periodic = (true, true),
    ),
    LSMSphericalShellDomain(;
        radius = 100.0,
        height = 3.5,
        nelements = (1, 10),
        npolynomial = 1,
    ),
]
orbital_data = Insolation.OrbitalData(joinpath(pkgdir(ClimaLSM), "artifacts"))
init_temp(z::FT, value::FT) where {FT} = FT(value)
for bucket_domain in bucket_domains
    @testset "Zero flux RHS" begin
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
            albedo,
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
        Y.bucket.T .= init_temp.(coords.subsurface.z, 280.0)
        Y.bucket.W .= 0.0 # no moisture
        Y.bucket.Ws .= 0.0 # no runoff
        Y.bucket.σS .= 0.0

        ode_function! = make_ode_function(model)
        set_initial_aux_state! = make_set_initial_aux_state(model)
        set_initial_aux_state!(p, Y, 0.0)
        dY = similar(Y)
        # init to nonzero numbers
        dY.bucket.T .= init_temp.(coords.subsurface.z, 1.0)
        dY.bucket.W .= 1.0
        dY.bucket.Ws .= 1.0
        dY.bucket.σS .= 0.0
        ode_function!(dY, Y, p, 0.0)
        @test mean(parent(dY.bucket.T)) < eps(FT)
        @test mean(parent(dY.bucket.W)) < eps(FT)
        @test mean(parent(dY.bucket.Ws)) < eps(FT)
        @test mean(parent(dY.bucket.σS)) < eps(FT)
    end


    @testset "Energy + Moisture Conservation" begin
        "Radiation"
        SW_d = (t) -> eltype(t)(10.0)
        LW_d = (t) -> eltype(t)(300.0)
        bucket_rad = PrescribedRadiativeFluxes(FT, SW_d, LW_d; orbital_data)
        "Atmos"
        precip = (t) -> eltype(t)(1e-6)
        T_atmos = (t) -> eltype(t)(298.0)
        u_atmos = (t) -> eltype(t)(4.0)
        q_atmos = (t) -> eltype(t)(0.01)
        h_atmos = FT(30)
        P_atmos = (t) -> eltype(t)(101325)
        bucket_atmos = PrescribedAtmosphere(
            precip,
            (t) -> eltype(t)(0.0),
            T_atmos,
            u_atmos,
            q_atmos,
            P_atmos,
            h_atmos,
        )
        Δt = FT(100.0)
        bucket_parameters = BucketModelParameters(
            κ_soil,
            ρc_soil,
            albedo,
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
        Y.bucket.T .= init_temp.(coords.subsurface.z, 280.0)
        Y.bucket.W .= 0.149
        Y.bucket.Ws .= 0.0
        Y.bucket.σS .= 0.0
        ode_function! = make_ode_function(model)
        t0 = 0.0
        set_initial_aux_state! = make_set_initial_aux_state(model)
        set_initial_aux_state!(p, Y, t0)
        dY = similar(Y)
        ode_function!(dY, Y, p, t0)
        F_water_sfc = precip(t0) .- p.bucket.evaporation
        F_sfc = -1 .* (p.bucket.turbulent_energy_flux .+ p.bucket.R_n)
        surface_space = model.domain.surface.space
        A_sfc = sum(ones(surface_space))
        # For the point space, we actually want the flux itself, since our variables are per unit area.
        # Divide by A_sfc
        if typeof(model.domain.surface.space) <: ClimaCore.Spaces.PointSpace
            F_sfc .= F_sfc ./ A_sfc
        end

        @test sum(F_water_sfc) ≈ sum(dY.bucket.W .+ dY.bucket.Ws)
        @test sum(F_sfc) ≈ sum(dY.bucket.T .* ρc_soil)

    end

end

using Test

using OrdinaryDiffEq: ODEProblem, solve, RK4, Euler
using DiffEqCallbacks
using DiffEqBase
using Statistics
using ClimaCore
import CLIMAParameters as CP

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM.Bucket:
    BucketModel,
    BucketModelParameters,
    PrescribedAtmosphere,
    PrescribedRadiativeFluxes,
    BulkAlbedo
using ClimaLSM.Domains:
    coordinates,
    LSMSingleColumnDomain,
    LSMMultiColumnDomain,
    LSMSphericalShellDomain
using ClimaLSM: initialize, make_update_aux, make_ode_function


FT = Float64

# Bucket model parameters
import ClimaLSM
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))
earth_param_set = create_lsm_parameters(FT)

α_soil = (coordinate_point) -> FT(0.2) # soil albedo, spatially constant
α_snow = FT(0.8) # snow albedo
albedo = BulkAlbedo{FT}(α_snow, α_soil)
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
init_temp(z::FT, value::FT) where {FT} = FT(value)
for bucket_domain in bucket_domains
    @testset "Zero flux RHS" begin
        # Radiation
        SW_d = (t) -> eltype(t)(0.0)
        LW_d = (t) -> eltype(t)(5.67e-8 * 280.0^4.0)
        bucket_rad = PrescribedRadiativeFluxes(FT, SW_d, LW_d)
        # Atmos
        precip = (t) -> eltype(t)(0) # no precipitation
        T_atmos = (t) -> eltype(t)(280.0)
        u_atmos = (t) -> eltype(t)(1.0)
        q_atmos = (t) -> eltype(t)(0.0) # no atmos water
        h_atmos = FT(1e-8)
        ρ_atmos = (t) -> eltype(t)(1.13)
        ρ_sfc = FT(1.15)
        bucket_atmos = PrescribedAtmosphere(
            precip,
            precip,
            T_atmos,
            u_atmos,
            q_atmos,
            ρ_atmos,
            h_atmos,
            ρ_sfc,
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
        bucket_rad = PrescribedRadiativeFluxes(FT, SW_d, LW_d)
        "Atmos"
        precip = (t) -> eltype(t)(1e-6)
        T_atmos = (t) -> eltype(t)(298.0)
        u_atmos = (t) -> eltype(t)(4.0)
        q_atmos = (t) -> eltype(t)(0.01)
        h_atmos = FT(30)
        ρ_atmos = (t) -> eltype(t)(1.13)
        ρ_sfc = FT(1.15)
        bucket_atmos = PrescribedAtmosphere(
            precip,
            (t) -> eltype(t)(0.0),
            T_atmos,
            u_atmos,
            q_atmos,
            ρ_atmos,
            h_atmos,
            ρ_sfc,
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
        tf = 10000.0
        prob = ODEProblem(ode_function!, Y, (t0, tf), p)

        # need a callback to get and store `p`
        saved_values = SavedValues(FT, ClimaCore.Fields.FieldVector)
        cb = SavingCallback(
            (u, t, integrator) -> copy(integrator.p),
            saved_values;
            saveat = 0:Δt:tf,
        )
        sol = solve(prob, Euler(); dt = Δt, saveat = 0:Δt:tf, callback = cb)

        # Sums integrate over area
        # ∫ variable dA
        W = [sum(sol.u[k].bucket.W) for k in 1:length(sol.t)]
        Ws = [sum(sol.u[k].bucket.Ws) for k in 1:length(sol.t)]
        R_n = [sum(saved_values.saveval[k].bucket.R_n) for k in 1:length(sol.t)]
        turbulent_energy_flux = [
            sum(saved_values.saveval[k].bucket.turbulent_energy_flux) for
            k in 1:length(sol.t)
        ]
        evaporation = [
            sum(saved_values.saveval[k].bucket.evaporation) for
            k in 1:length(sol.t)
        ]

        F_sfc = turbulent_energy_flux .+ R_n
        surface_space = model.domain.surface.space
        A_sfc = sum(ones(surface_space))
        F_water_sfc = precip.(sol.t) * A_sfc .- evaporation
        net_water = W .+ Ws
        # dY(t+dt) = Y(t+dt) - W(t) = RHS|_t * dt for Euler stepping
        # First element of `p` not filled in, so start at 2.
        # If we call update aux before the simulation, it should be.

        # Look at fractional error per step

        # For the point space, we actually want the flux itself, since our energy is per unit area.
        # Divide by A_sfc
        if typeof(model.domain.surface.space) <: ClimaCore.Spaces.PointSpace
            F_sfc .= F_sfc ./ A_sfc
        end

        dE_land = ρc_soil * (sum(sol.u[end].bucket.T) - sum(sol.u[2].bucket.T))
        dE_expected = -sum(F_sfc[2:(end - 1)]) * Δt
        @test abs(dE_land - dE_expected) ./
              (ρc_soil * sum(sol.u[2].bucket.T) ./ length(sol.t)) < 1e-13
        dW_land = net_water[end] - net_water[2]
        dW_expected = sum(F_water_sfc[2:(end - 1)]) * Δt

        @test abs(dW_land - dW_expected) ./ net_water[2] ./ length(sol.t) <
              1e-13

    end

end

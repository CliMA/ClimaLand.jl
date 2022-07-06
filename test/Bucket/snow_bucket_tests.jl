using Test

using OrdinaryDiffEq: ODEProblem, solve, Euler
using DiffEqCallbacks
using Statistics
using ClimaCore

if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM.Bucket:
    BucketModel,
    BucketModelParameters,
    PrescribedAtmosphere,
    PrescribedRadiativeFluxes,
    BulkAlbedo,
    partition_surface_fluxes
using ClimaLSM.Domains: coordinates, LSMSingleColumnDomain, LSMMultiColumnDomain
using ClimaLSM: initialize, make_update_aux, make_ode_function

FT = Float64

# Bucket model parameters
import ClimaLSM
import ClimaLSM.Parameters as LSMP
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
]
init_temp(z::FT, value::FT) where {FT} = FT(value)
for bucket_domain in bucket_domains

    @testset "Conservation of water and energy" begin
        "Radiation"
        SW_d = (t) -> eltype(t)(20.0)
        LW_d = (t) -> eltype(t)(20.0)
        bucket_rad = PrescribedRadiativeFluxes(FT, SW_d, LW_d)
        "Atmos"
        precip = (t) -> eltype(t)(0) # no precipitation
        T_atmos = (t) -> eltype(t)(280.0)
        u_atmos = (t) -> eltype(t)(10.0)
        q_atmos = (t) -> eltype(t)(0.03)
        h_atmos = FT(3)
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
        τc = FT(10.0)
        bucket_parameters = BucketModelParameters(
            κ_soil,
            ρc_soil,
            albedo,
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

        # run for some time
        Y, p, coords = initialize(model)
        Y.bucket.T .= init_temp.(coords.subsurface.z, 274.0)
        Y.bucket.W .= 0.0 # no moisture
        Y.bucket.Ws .= 0.0 # no runoff
        Y.bucket.σS .= 0.5

        ode_function! = make_ode_function(model)
        update_aux! = make_update_aux(model)
        update_aux!(p, Y, 0.0)
        saved_values = SavedValues(FT, ClimaCore.Fields.FieldVector)
        cb = SavingCallback(
            (u, t, integrator) -> copy(integrator.p),
            saved_values;
            saveat = 0:Δt:100.0,
        )
        prob = ODEProblem(ode_function!, Y, (0.0, 100.0), p)
        sol = solve(
            prob,
            Euler();
            dt = Δt,
            reltol = 1e-6,
            abstol = 1e-6,
            callback = cb,
        )

        _LH_f0 = LSMP.LH_f0(model.parameters.earth_param_set)
        _T_freeze = LSMP.T_freeze(model.parameters.earth_param_set)
        snow_cover_fraction(σS) = σS > eps(FT) ? FT(1.0) : FT(0.0)

        # Extract water content on land and fluxes
        W = [sum(sol.u[k].bucket.W) for k in 1:length(sol.t)]
        σS = [sum(sol.u[k].bucket.σS) for k in 1:length(sol.t)]
        Ws = [sum(sol.u[k].bucket.Ws) for k in 1:length(sol.t)]

        R_n = [sum(saved_values.saveval[k].bucket.R_n) for k in 1:length(sol.t)] # net Radiation, ∫dA
        turbulent_energy_flux = [
            sum(saved_values.saveval[k].bucket.turbulent_energy_flux) for
            k in 1:length(sol.t)
        ] # ∫ dA
        evaporation = [
            sum(saved_values.saveval[k].bucket.evaporation) for
            k in 1:length(sol.t)
        ] # ∫ dA

        partitioned_fluxes = [
            partition_surface_fluxes.(
                sol.u[k].bucket.σS,
                saved_values.saveval[k].bucket.T_sfc,
                model.parameters.τc,
                snow_cover_fraction.(sol.u[k].bucket.σS),
                saved_values.saveval[k].bucket.evaporation,
                saved_values.saveval[k].bucket.turbulent_energy_flux .+
                saved_values.saveval[k].bucket.R_n,
                _LH_f0,
                _T_freeze,
            ) for k in 1:length(sol.t)
        ]
        F_melt = [
            sum(partitioned_fluxes[k].F_melt) for
            k in 1:length(partitioned_fluxes)
        ]
        F_into_snow = [
            sum(partitioned_fluxes[k].F_into_snow) for
            k in 1:length(partitioned_fluxes)
        ]
        G = [sum(partitioned_fluxes[k].G) for k in 1:length(partitioned_fluxes)]

        F_sfc = turbulent_energy_flux .+ R_n
        surface_space = model.domain.surface.space
        A_sfc = sum(ones(surface_space))
        F_water_sfc = precip.(sol.t) * A_sfc .- evaporation

        if typeof(model.domain.surface.space) <: ClimaCore.Spaces.PointSpace
            G .= G ./ A_sfc
            F_sfc .= F_sfc ./ A_sfc
            F_into_snow .= F_into_snow ./ A_sfc
            σS .= σS ./ A_sfc
            W .= W ./ A_sfc
            Ws .= Ws ./ A_sfc
            evaporation .= evaporation ./ A_sfc
        end

        # compute total energy and water contents
        e_soil = [sum(sol.u[k].bucket.T) for k in 1:length(sol.t)] .* ρc_soil
        Isnow = -_LH_f0 * σS
        WL = W .+ σS .+ Ws
        IL = Isnow .+ e_soil

        # check time derivatives of components and total:
        dIsnowdt = (Isnow[2:end] - Isnow[1:(end - 1)]) / Δt
        # dIsnow/dt = -F_into_snow
        Δsnow = dIsnowdt .+ F_into_snow[1:(end - 1)]
        @test maximum(abs.(Δsnow[2:end])) ./ Isnow[2] * Δt < 10.0 * eps(FT)
        # de_soil/dt = -G
        dedt = (e_soil[2:end] - e_soil[1:(end - 1)]) / Δt
        Δe = dedt .+ G[1:(end - 1)]
        @test maximum(abs.(Δe[2:end])) ./ e_soil[2] * Δt < 10.0 * eps(FT)
        # dW/dt = -E -> dWdt + E = 0
        dWdt = (WL[2:end] - WL[1:(end - 1)]) / Δt
        @test maximum(abs.(dWdt[2:end] .+ evaporation[2:(end - 1)])) * Δt /
              WL[2] < 10.0 * eps(FT)
        # dIL/dt = - (F_sfc .- 0.0)
        dILdt = (IL[2:end] - IL[1:(end - 1)]) / Δt
        Δ = dILdt .+ F_sfc[1:(end - 1)]
        @test maximum(abs.(Δ[2:end])) ./ IL[2] * Δt < 10.0 * eps(FT)

    end
end

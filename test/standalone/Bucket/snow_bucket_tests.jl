using Test

using Statistics
using Insolation
using Dates
using ClimaCore
using ClimaLSM.Bucket:
    BucketModel,
    BucketModelParameters,
    BulkAlbedoFunction,
    partition_surface_fluxes
using ClimaLSM.Domains: coordinates, Column, HybridBox, SphericalShell
using ClimaLSM:
    initialize,
    make_update_aux,
    make_exp_tendency,
    make_set_initial_aux_state,
    PrescribedAtmosphere,
    PrescribedRadiativeFluxes

FT = Float64

# Bucket model parameters
import ClimaLSM
import ClimaLSM.Parameters as LSMP
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
    Column(; zlim = (-100.0, 0.0), nelements = 10),
    HybridBox(;
        xlim = (-1.0, 0.0),
        ylim = (-1.0, 0.0),
        zlim = (-100.0, 0.0),
        nelements = (2, 2, 10),
        npolynomial = 1,
        periodic = (true, true),
    ),
    SphericalShell(;
        radius = 100.0,
        depth = 3.5,
        nelements = (1, 10),
        npolynomial = 1,
    ),
]
init_temp(z::FT, value::FT) where {FT} = FT(value)
orbital_data = Insolation.OrbitalData()
for i in 1:3
    @testset "Conservation of water and energy I (snow present)" begin
        "Radiation"
        ref_time = DateTime(2005)
        SW_d = (t) -> eltype(t)(20.0)
        LW_d = (t) -> eltype(t)(20.0)
        bucket_rad =
            PrescribedRadiativeFluxes(FT, SW_d, LW_d, ref_time; orbital_data)
        "Atmos"
        liquid_precip = (t) -> eltype(t)(1e-8) # precipitation
        snow_precip = (t) -> eltype(t)(1e-7) # precipitation

        T_atmos = (t) -> eltype(t)(280.0)
        u_atmos = (t) -> eltype(t)(10.0)
        q_atmos = (t) -> eltype(t)(0.03)
        h_atmos = FT(3)
        P_atmos = (t) -> eltype(t)(101325) # Pa
        bucket_atmos = PrescribedAtmosphere(
            liquid_precip,
            snow_precip,
            T_atmos,
            u_atmos,
            q_atmos,
            P_atmos,
            ref_time,
            h_atmos,
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
            domain = bucket_domains[i],
            atmosphere = bucket_atmos,
            radiation = bucket_rad,
        )

        # run for some time
        Y, p, coords = initialize(model)
        Y.bucket.T .= init_temp.(coords.subsurface.z, 274.0)
        Y.bucket.W .= 0.0 # no moisture
        Y.bucket.Ws .= 0.0 # no runoff
        Y.bucket.σS .= 0.5
        t0 = FT(0.0)
        dY = similar(Y)

        exp_tendency! = make_exp_tendency(model)
        set_initial_aux_state! = make_set_initial_aux_state(model)
        set_initial_aux_state!(p, Y, t0)
        exp_tendency!(dY, Y, p, t0)
        _LH_f0 = LSMP.LH_f0(model.parameters.earth_param_set)
        _ρ_liq = LSMP.ρ_cloud_liq(model.parameters.earth_param_set)
        _ρLH_f0 = _ρ_liq * _LH_f0 # Latent heat per unit volume
        _T_freeze = LSMP.T_freeze(model.parameters.earth_param_set)
        snow_cover_fraction(σS) = σS > eps(FT) ? FT(1.0) : FT(0.0)

        partitioned_fluxes =
            partition_surface_fluxes.(
                Y.bucket.σS,
                p.bucket.T_sfc,
                model.parameters.τc,
                snow_cover_fraction.(Y.bucket.σS),
                p.bucket.evaporation,
                p.bucket.turbulent_energy_flux .+ p.bucket.R_n,
                _ρLH_f0,
                _T_freeze,
            )
        F_melt = partitioned_fluxes.F_melt
        F_into_snow =
            partitioned_fluxes.F_into_snow .+ _ρLH_f0 .* snow_precip(t0)
        G = partitioned_fluxes.G
        F_sfc =
            p.bucket.turbulent_energy_flux .+ p.bucket.R_n .+
            _ρLH_f0 .* snow_precip(t0)
        F_water_sfc =
            liquid_precip(t0) + snow_precip(t0) .- p.bucket.evaporation


        if i == 1
            A_point = sum(ones(bucket_domains[i].space.surface))
        else
            A_point = 1
        end

        dIsnow = -_ρLH_f0 .* dY.bucket.σS
        @test sum(dIsnow) / A_point ≈ sum(-1 .* F_into_snow) / A_point

        de_soil = dY.bucket.T .* ρc_soil
        @test sum(de_soil) ≈ sum(-1 .* G) / A_point

        dWL = dY.bucket.W .+ dY.bucket.Ws .+ dY.bucket.σS
        @test sum(dWL) / A_point ≈ sum(F_water_sfc) / A_point

        dIL = sum(dIsnow) / A_point .+ sum(de_soil)
        @test dIL ≈ sum(-1 .* F_sfc) / A_point
    end
end

for i in 1:3
    @testset "Conservation of water and energy II (no snow to start)" begin
        "Radiation"
        ref_time = DateTime(2005)
        SW_d = (t) -> eltype(t)(20.0)
        LW_d = (t) -> eltype(t)(20.0)
        bucket_rad =
            PrescribedRadiativeFluxes(FT, SW_d, LW_d, ref_time; orbital_data)
        "Atmos"
        liquid_precip = (t) -> eltype(t)(1e-8) # precipitation
        snow_precip = (t) -> eltype(t)(1e-7) # precipitation

        T_atmos = (t) -> eltype(t)(280.0)
        u_atmos = (t) -> eltype(t)(10.0)
        q_atmos = (t) -> eltype(t)(0.03)
        h_atmos = FT(3)
        P_atmos = (t) -> eltype(t)(101325) # Pa
        bucket_atmos = PrescribedAtmosphere(
            liquid_precip,
            snow_precip,
            T_atmos,
            u_atmos,
            q_atmos,
            P_atmos,
            ref_time,
            h_atmos,
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
            domain = bucket_domains[i],
            atmosphere = bucket_atmos,
            radiation = bucket_rad,
        )

        # run for some time
        Y, p, coords = initialize(model)
        Y.bucket.T .= init_temp.(coords.subsurface.z, 274.0)
        Y.bucket.W .= 0.0 # no moisture
        Y.bucket.Ws .= 0.0 # no runoff
        Y.bucket.σS .= 0.0
        t0 = FT(0.0)
        dY = similar(Y)

        exp_tendency! = make_exp_tendency(model)
        set_initial_aux_state! = make_set_initial_aux_state(model)
        set_initial_aux_state!(p, Y, t0)
        exp_tendency!(dY, Y, p, t0)
        _LH_f0 = LSMP.LH_f0(model.parameters.earth_param_set)
        _ρ_liq = LSMP.ρ_cloud_liq(model.parameters.earth_param_set)
        _ρLH_f0 = _ρ_liq * _LH_f0 # Latent heat per unit volume
        _T_freeze = LSMP.T_freeze(model.parameters.earth_param_set)
        snow_cover_fraction(σS) = σS > eps(FT) ? FT(1.0) : FT(0.0)

        partitioned_fluxes =
            partition_surface_fluxes.(
                Y.bucket.σS,
                p.bucket.T_sfc,
                model.parameters.τc,
                snow_cover_fraction.(Y.bucket.σS),
                p.bucket.evaporation,
                p.bucket.turbulent_energy_flux .+ p.bucket.R_n,
                _ρLH_f0,
                _T_freeze,
            )
        F_melt = partitioned_fluxes.F_melt
        F_into_snow =
            partitioned_fluxes.F_into_snow .+ _ρLH_f0 .* snow_precip(t0)
        G = partitioned_fluxes.G
        F_sfc =
            p.bucket.turbulent_energy_flux .+ p.bucket.R_n .+
            _ρLH_f0 .* snow_precip(t0)
        F_water_sfc =
            liquid_precip(t0) + snow_precip(t0) .- p.bucket.evaporation

        if i == 1
            A_point = sum(ones(bucket_domains[i].space.surface))
        else
            A_point = 1
        end

        dIsnow = -_ρLH_f0 .* dY.bucket.σS
        @test sum(dIsnow) / A_point ≈ sum(-1 .* F_into_snow) / A_point

        de_soil = dY.bucket.T .* ρc_soil
        @test sum(de_soil) ≈ sum(-1 .* G) / A_point

        dWL = dY.bucket.W .+ dY.bucket.Ws .+ dY.bucket.σS
        @test sum(dWL) / A_point ≈ sum(F_water_sfc) / A_point

        dIL = sum(dIsnow) / A_point .+ sum(de_soil)
        @test dIL ≈ sum(-1 .* F_sfc) / A_point
    end
end



@testset "Conservation of water and energy - nonuniform evaporation" begin
    i = 3
    "Radiation"
    ref_time = DateTime(2005)
    SW_d = (t) -> eltype(t)(20.0)
    LW_d = (t) -> eltype(t)(20.0)
    bucket_rad =
        PrescribedRadiativeFluxes(FT, SW_d, LW_d, ref_time; orbital_data)
    "Atmos"
    liquid_precip = (t) -> eltype(t)(1e-8) # precipitation
    snow_precip = (t) -> eltype(t)(1e-7) # precipitation

    T_atmos = (t) -> eltype(t)(280.0)
    u_atmos = (t) -> eltype(t)(10.0)
    q_atmos = (t) -> eltype(t)(0.03)
    h_atmos = FT(3)
    P_atmos = (t) -> eltype(t)(101325) # Pa
    bucket_atmos = PrescribedAtmosphere(
        liquid_precip,
        snow_precip,
        T_atmos,
        u_atmos,
        q_atmos,
        P_atmos,
        ref_time,
        h_atmos,
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
        domain = bucket_domains[i],
        atmosphere = bucket_atmos,
        radiation = bucket_rad,
    )

    # run for some time
    Y, p, coords = initialize(model)
    Y.bucket.T .= init_temp.(coords.subsurface.z, 274.0)
    Y.bucket.W .= 0.0 # no moisture
    Y.bucket.Ws .= 0.0 # no runoff
    Y.bucket.σS .= 0.0
    t0 = FT(0.0)
    dY = similar(Y)

    compute_exp_tendency! = ClimaLSM.make_compute_exp_tendency(model)
    set_initial_aux_state! = make_set_initial_aux_state(model)
    set_initial_aux_state!(p, Y, t0)
    random = zeros(bucket_domains[i].space.surface)
    parent(random) .= rand(FT, size(parent(random)))
    p.bucket.evaporation .= random .* 1e-7
    compute_exp_tendency!(dY, Y, p, t0)
    _LH_f0 = LSMP.LH_f0(model.parameters.earth_param_set)
    _ρ_liq = LSMP.ρ_cloud_liq(model.parameters.earth_param_set)
    _ρLH_f0 = _ρ_liq * _LH_f0 # Latent heat per unit volume
    _T_freeze = LSMP.T_freeze(model.parameters.earth_param_set)
    snow_cover_fraction(σS) = σS > eps(FT) ? FT(1.0) : FT(0.0)

    partitioned_fluxes =
        partition_surface_fluxes.(
            Y.bucket.σS,
            p.bucket.T_sfc,
            model.parameters.τc,
            snow_cover_fraction.(Y.bucket.σS),
            p.bucket.evaporation,
            p.bucket.turbulent_energy_flux .+ p.bucket.R_n,
            _ρLH_f0,
            _T_freeze,
        )
    F_melt = partitioned_fluxes.F_melt
    F_into_snow = partitioned_fluxes.F_into_snow .+ _ρLH_f0 .* snow_precip(t0)
    G = partitioned_fluxes.G
    F_sfc =
        p.bucket.turbulent_energy_flux .+ p.bucket.R_n .+
        _ρLH_f0 .* snow_precip(t0)
    F_water_sfc = liquid_precip(t0) + snow_precip(t0) .- p.bucket.evaporation

    dIsnow = -_ρLH_f0 .* dY.bucket.σS
    @test sum(dIsnow) ≈ sum(-1 .* F_into_snow)

    de_soil = dY.bucket.T .* ρc_soil
    @test sum(de_soil) ≈ sum(-1 .* G)

    dWL = dY.bucket.W .+ dY.bucket.Ws .+ dY.bucket.σS
    @test sum(dWL) ≈ sum(F_water_sfc)

    dIL = sum(dIsnow) .+ sum(de_soil)
    @test dIL ≈ sum(-1 .* F_sfc)
end

using Test

using Thermodynamics
using SurfaceFluxes
using StaticArrays
using ClimaLSM.Snow
using ClimaLSM.Domains
import ClimaLSM
import ClimaLSM.Parameters as LSMP
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))

@testset "Snow Model" begin
    FT = Float32
    param_set = create_lsm_parameters(FT)
    Δt = FT(180.0)
    parameters = SnowParameters{FT}(Δt; earth_param_set = param_set)
    domain = Point(; z_sfc = FT(0))
    "Radiation"
    SW_d = (t) -> eltype(t)(20.0)
    LW_d = (t) -> eltype(t)(20.0)
    rad = ClimaLSM.PrescribedRadiativeFluxes(FT, SW_d, LW_d)
    "Atmos"
    precip = (t) -> eltype(t)(0) # no precipitation
    T_atmos = (t) -> eltype(t)(290.0)
    u_atmos = (t) -> eltype(t)(10.0)
    q_atmos = (t) -> eltype(t)(0.03)
    h_atmos = FT(3)
    ρ_atmos = (t) -> eltype(t)(1.13)
    atmos = ClimaLSM.PrescribedAtmosphere(
        precip,
        precip,
        T_atmos,
        u_atmos,
        q_atmos,
        ρ_atmos,
        h_atmos,
    )
    model = ClimaLSM.Snow.SnowModel(
        parameters = parameters,
        domain = domain,
        atmosphere = atmos,
        radiation = rad,
    )
    Y, p, coords = ClimaLSM.initialize(model)
    @test (Y.snow |> propertynames) == (:S, :U)
    @test (p.snow |> propertynames) == (
        :q_l,
        :T,
        :T_sfc,
        :ρ_sfc,
        :q_sfc,
        :evaporation,
        :turbulent_energy_flux,
        :R_n,
        :energy_runoff,
        :water_runoff,
        :total_energy_flux,
        :total_water_flux,
        :snow_energy_flux,
        :snow_water_flux,
    )

    Y.snow.S .= FT(0.1)
    Y.snow.U .=
        ClimaLSM.Snow.energy_from_temperature_and_swe.(
            Y.snow.S,
            FT(273.0),
            Ref(model.parameters),
        )
    update_aux! = ClimaLSM.make_update_aux(model)
    t0 = FT(0.0)
    update_aux!(p, Y, t0)
    (; α_snow, ϵ_snow) = model.parameters
    _σ = LSMP.Stefan(model.parameters.earth_param_set)

    # Check if aux update occurred correctly
    @test p.snow.R_n == @.(
        -(1 - α_snow) * SW_d(t0) - ϵ_snow * (LW_d(t0) - _σ * p.snow.T_sfc^4)
    )
    @test p.snow.R_n == ClimaLSM.net_radiation(model.radiation, model, Y, p, t0)
    @test p.snow.q_l ==
          snow_liquid_mass_fraction.(Y.snow.U, Y.snow.S, Ref(model.parameters))
    @test p.snow.T_sfc ==
          snow_surface_temperature.(
        snow_bulk_temperature.(
            Y.snow.U,
            Y.snow.S,
            p.snow.q_l,
            Ref(model.parameters),
        )
    )

    ρ_sfc =
        ClimaLSM.surface_air_density(model.atmos, model, Y, p, t0, p.snow.T_sfc)
    thermo_params =
        LSMP.thermodynamic_parameters(model.parameters.earth_param_set)
    q_sfc =
        Thermodynamics.q_vap_saturation_generic.(
            Ref(thermo_params),
            p.snow.T_sfc,
            ρ_sfc,
            Ref(Thermodynamics.Ice()),
        )
    turbulent_fluxes = ClimaLSM.surface_fluxes(model.atmos, model, Y, p, t0)
    @test turbulent_fluxes.turbulent_energy_flux == p.snow.turbulent_energy_flux
    @test turbulent_fluxes.evaporation == p.snow.evaporation

    # Now compute tendencies and make sure they operate correctly.
    dY = similar(Y)
    exp_tendency! = ClimaLSM.make_compute_exp_tendency(model)
    exp_tendency!(dY, Y, p, t0)
    _ρ_liq = LSMP.ρ_cloud_liq(model.parameters.earth_param_set)
    net_water_fluxes = @.(
        -model.atmos.liquid_precip(t0) - model.atmos.snow_precip(t0) -
        p.snow.evaporation + p.snow.water_runoff
    )
    @test dY.snow.S == net_water_fluxes
    @test dY.snow.U ==
          @.(-p.snow.turbulent_energy_flux - p.snow.R_n + p.snow.energy_runoff)


    # Now try a step where the snow will melt
    Y.snow.S .= FT(0.1)
    Y.snow.U .= FT(-2e6)
    update_aux!(p, Y, t0)
    dY = similar(Y)
    exp_tendency!(dY, Y, p, t0)
    @test dY.snow.S == @.(-Y.snow.S / Δt)
    @test all(
        ClimaLSM.Snow.clip_dSdt.([0.1, 1.0, 2.0], [-1.0, -1.0, -1.0], 1.0) .==
        [-0.1, -1.0, -1.0],
    )
end

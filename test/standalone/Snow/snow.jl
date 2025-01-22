using Test
import ClimaComms
ClimaComms.@import_required_backends
import ClimaParams as CP
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
using Thermodynamics
using SurfaceFluxes
using StaticArrays
using Dates
import ClimaLand
using ClimaLand.Snow
using ClimaLand.Domains
import ClimaLand.Parameters as LP

@testset "Snow Model" begin
    FT = Float32
    earth_param_set = LP.LandParameters(FT)

    start_date = DateTime(2005)
    param_set = LP.LandParameters(FT)
    Δt = FT(180.0)
    parameters = SnowParameters{FT}(Δt; earth_param_set = param_set)
    domain = Point(; z_sfc = FT(0))
    "Radiation"
    SW_d = TimeVaryingInput((t) -> eltype(t)(20.0))
    LW_d = TimeVaryingInput((t) -> eltype(t)(20.0))
    rad = ClimaLand.PrescribedRadiativeFluxes(FT, SW_d, LW_d, start_date)
    "Atmos"
    precip = TimeVaryingInput((t) -> eltype(t)(0)) # no precipitation
    T_atmos = TimeVaryingInput((t) -> eltype(t)(290.0))
    u_atmos = TimeVaryingInput((t) -> eltype(t)(10.0))
    q_atmos = TimeVaryingInput((t) -> eltype(t)(0.03))
    h_atmos = FT(3)
    P_atmos = TimeVaryingInput((t) -> eltype(t)(101325))
    atmos = ClimaLand.PrescribedAtmosphere(
        precip,
        precip,
        T_atmos,
        u_atmos,
        q_atmos,
        P_atmos,
        start_date,
        h_atmos,
        earth_param_set,
    )
    model = ClimaLand.Snow.SnowModel(
        parameters = parameters,
        domain = domain,
        boundary_conditions = ClimaLand.Snow.AtmosDrivenSnowBC(atmos, rad),
    )
    drivers = ClimaLand.get_drivers(model)
    @test drivers == (atmos, rad)
    Y, p, coords = ClimaLand.initialize(model)
    @test (Y.snow |> propertynames) == (:S, :U)
    @test (p.snow |> propertynames) == (
        :q_sfc,
        :q_l,
        :κ,
        :T,
        :T_sfc,
        :ρ_snow,
        :turbulent_fluxes,
        :R_n,
        :energy_runoff,
        :water_runoff,
        :total_energy_flux,
        :total_water_flux,
        :applied_energy_flux,
        :applied_water_flux,
        :snow_cover_fraction,
    )

    Y.snow.S .= FT(0.1)
    Y.snow.U .=
        ClimaLand.Snow.energy_from_T_and_swe.(
            Y.snow.S,
            FT(273.0),
            Ref(model.parameters),
        )
    set_initial_cache! = ClimaLand.make_set_initial_cache(model)
    t0 = FT(0.0)
    set_initial_cache!(p, Y, t0)
    (; α_snow, ϵ_snow) = model.parameters
    _σ = LP.Stefan(model.parameters.earth_param_set)

    # Check if aux update occurred correctly
    @test p.snow.R_n ==
          @. (-(1 - α_snow) * 20.0f0 - ϵ_snow * (20.0f0 - _σ * p.snow.T_sfc^4))
    R_n_copy = copy(p.snow.R_n)
    ClimaLand.net_radiation!(
        R_n_copy,
        model.boundary_conditions.radiation,
        model,
        Y,
        p,
        t0,
    )
    @test p.snow.R_n == R_n_copy
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

    ρ_sfc = ClimaLand.surface_air_density(
        model.boundary_conditions.atmos,
        model,
        Y,
        p,
        t0,
        p.snow.T_sfc,
    )
    thermo_params =
        LP.thermodynamic_parameters(model.parameters.earth_param_set)
    q_sfc =
        Thermodynamics.q_vap_saturation_generic.(
            Ref(thermo_params),
            p.snow.T_sfc,
            ρ_sfc,
            Ref(Thermodynamics.Ice()),
        )
    @test p.snow.q_sfc ≈ q_sfc
    turb_fluxes_copy = copy(p.snow.turbulent_fluxes)
    ClimaLand.turbulent_fluxes!(
        turb_fluxes_copy,
        model.boundary_conditions.atmos,
        model,
        Y,
        p,
        t0,
    )
    @test turb_fluxes_copy.shf == p.snow.turbulent_fluxes.shf
    @test turb_fluxes_copy.lhf == p.snow.turbulent_fluxes.lhf
    @test turb_fluxes_copy.vapor_flux == p.snow.turbulent_fluxes.vapor_flux
    old_ρ = deepcopy(p.snow.ρ_snow)
    Snow.update_density!(model.parameters.density, model.parameters, Y, p)
    @test p.snow.ρ_snow == old_ρ
    old_z = similar(Y.snow.S)
    old_z .= FT(0.5)
    z = snow_depth(model.parameters.density, Y, p, parameters)
    @test z == old_z

    # Now compute tendencies and make sure they operate correctly.
    dY = similar(Y)
    exp_tendency! = ClimaLand.make_compute_exp_tendency(model)
    exp_tendency!(dY, Y, p, t0)
    _ρ_liq = LP.ρ_cloud_liq(model.parameters.earth_param_set)
    net_water_fluxes = @.(
        -0.0f0 - 0.0f0 - p.snow.turbulent_fluxes.vapor_flux +
        p.snow.water_runoff
    )
    @test dY.snow.S == net_water_fluxes
    @test dY.snow.U == @.(
        -p.snow.turbulent_fluxes.shf - p.snow.turbulent_fluxes.lhf -
        p.snow.R_n + p.snow.energy_runoff
    )
    @test isnothing(
        Snow.update_density_prog!(model.parameters.density, model, dY, Y, p),
    )


    # Now try a step where the snow will melt
    Y.snow.S .= FT(0.1)
    Y.snow.U .= FT(-2e6)
    set_initial_cache!(p, Y, t0)
    dY = similar(Y)
    exp_tendency!(dY, Y, p, t0)
    @test dY.snow.S == @.(-Y.snow.S / Δt)
    @test all(
        ClimaLand.Snow.clip_total_snow_water_flux.(
            [0.1, 1.0, 2.0],
            [1.0, 1.0, 1.0],
            1.0,
        ) .== [0.1, 1.0, 1.0],
    )
end

using Test
using Dates
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
import ClimaParams as CP
using ClimaLand
using ClimaLand.Snow
using ClimaLand.Domains: Column, HybridBox

import ClimaLand
import ClimaLand.Parameters as LP

for FT in (Float32, Float64)
    default_params_filepath =
        joinpath(pkgdir(ClimaLand), "toml", "default_parameters.toml")
    toml_dict = LP.create_toml_dict(FT, default_params_filepath)
    cmax = FT(0)
    cmin = FT(-2)
    nelems = 10
    col = Column(; zlim = (cmin, cmax), nelements = nelems)
    box = ClimaLand.Domains.HybridBox(;
        xlim = (cmin, cmax),
        ylim = (cmin, cmax),
        zlim = (cmin, cmax),
        nelements = (nelems, nelems, nelems),
    )

    domains = [col, box]

    t0 = 0.0

    earth_param_set = LP.LandParameters(toml_dict)

    start_date = DateTime(2005)
    Δt = FT(180.0)
    parameters = SnowParameters(toml_dict, Δt)
    "Radiation"
    SW_d = TimeVaryingInput((t) -> eltype(t)(20.0))
    LW_d = TimeVaryingInput((t) -> eltype(t)(20.0))
    rad = ClimaLand.PrescribedRadiativeFluxes(FT, SW_d, LW_d, start_date)
    "Atmos"
    precip = TimeVaryingInput((t) -> eltype(t)(0)) # no precipitation
    T_atmos = TimeVaryingInput((t) -> eltype(t)(290.0))
    u_atmos = TimeVaryingInput((t) -> eltype(t)(10.0))
    q_atmos = TimeVaryingInput((t) -> eltype(t)(0.0003))
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

    @testset "Snow model total energy and water, FT = $FT" begin
        for domain in domains
            snow = Snow.SnowModel(
                parameters = parameters,
                domain = ClimaLand.Domains.obtain_surface_domain(domain),
                boundary_conditions = Snow.AtmosDrivenSnowBC(atmos, rad),
            )

            Y, p, cds = initialize(snow)
            S0 = FT(0.5)
            S_l0 = FT(0.3)
            Temp0 = FT(270)
            U0 = Snow.energy_from_T_and_swe.(S0, Temp0, snow.parameters)

            Y.snow.S .= S0
            Y.snow.S_l .= S_l0
            Y.snow.U .= U0

            set_initial_cache! = make_set_initial_cache(snow)
            set_initial_cache!(p, Y, t0)
            total_water = ClimaCore.Fields.zeros(domain.space.surface)
            ClimaLand.total_liq_water_vol_per_area!(total_water, snow, Y, p, t0)
            @test all(parent(total_water) .≈ S0)
            total_energy = ClimaCore.Fields.zeros(domain.space.surface)
            ClimaLand.total_energy_per_area!(total_energy, snow, Y, p, t0)
            @test all(parent(total_energy) .≈ U0)
        end
    end
end

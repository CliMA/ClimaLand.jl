using Test
import ClimaComms
ClimaComms.@import_required_backends
using Statistics
using ClimaCore
import ClimaParams as CP
using ClimaLand
using ClimaLand.Soil
using ClimaLand.Domains: Column, HybridBox

import ClimaLand
import ClimaLand.Parameters as LP
struct FakeSource{FT} <: AbstractSoilSource{FT} end

import ClimaLand.source!
# These allocate, but we dont mind here.
function ClimaLand.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::FakeSource{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    model::RichardsModel{FT},
) where {FT}
    source = ClimaLand.Fields.zeros(model.domain.space.subsurface) .- FT(1e-5)
    @. dY.soil.ϑ_l += source
    tmp = ClimaCore.Fields.zeros(model.domain.space.surface)
    ClimaCore.Operators.column_integral_definite!(tmp, source)
    @. p.soil.∫Swdz += tmp
end

function ClimaLand.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::FakeSource{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    model::EnergyHydrology{FT},
) where {FT}
    earth_param_set = model.parameters.earth_param_set
    _ρ_l = LP.ρ_cloud_liq(earth_param_set)
    _ρ_i = LP.ρ_cloud_ice(earth_param_set)
    source = ClimaLand.Fields.zeros(model.domain.space.subsurface) .- FT(1e-5)
    @. dY.soil.ϑ_l += source
    @. dY.soil.θ_i += source * (_ρ_l / _ρ_i)
    tmp = ClimaCore.Fields.zeros(model.domain.space.surface)
    ClimaCore.Operators.column_integral_definite!(tmp, source)
    @. p.soil.∫Swdz += tmp
    @. p.soil.∫Sedz = 0
end


for FT in (Float32, Float64)
    cmax = FT(0)
    cmin = FT(-2)
    nelems = 10
    col = Column(; zlim = (cmin, cmax), nelements = nelems)
    box = ClimaLand.Domains.HybridBox(;
        xlim = (cmin, cmax),
        ylim = (cmin, cmax),
        zlim = (cmin, cmax),
        nelements = (nelems, nelems, nelems),
        npolynomial = 0,
    )

    domains = [col, box]
    ν = FT(0.495)
    K_sat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(2.0)
    vg_α = FT(2.6) # inverse meters
    hcm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
    θ_r = FT(0)
    ν_ss_om = FT(0.0)
    ν_ss_quartz = FT(1.0)
    ν_ss_gravel = FT(0.0)
    top_water_flux_bc = WaterFluxBC((p, t) -> 1.0)
    top_heat_flux_bc = HeatFluxBC((p, t) -> 1.0)
    bot_water_flux_bc = WaterFluxBC((p, t) -> -1.0)
    bot_heat_flux_bc = HeatFluxBC((p, t) -> -1.0)
    t0 = 0.0
    @testset "Richards total water, FT = $FT" begin
        sources = (FakeSource{FT}(),)
        params = Soil.RichardsParameters(;
            ν = ν,
            hydrology_cm = hcm,
            K_sat = K_sat,
            S_s = S_s,
            θ_r = θ_r,
        )
        boundary_fluxes =
            (; top = top_water_flux_bc, bottom = bot_water_flux_bc)
        for domain in domains
            soil = Soil.RichardsModel{FT}(;
                parameters = params,
                domain = domain,
                boundary_conditions = boundary_fluxes,
                sources = sources,
            )

            Y, p, cds = initialize(soil)
            Y.soil.ϑ_l .= ν / 2
            set_initial_cache! = make_set_initial_cache(soil)
            set_initial_cache!(p, Y, t0)
            total_water = ClimaCore.Fields.zeros(soil.domain.space.surface)
            ClimaLand.total_liq_water_vol_per_area!(total_water, soil, Y, p, t0)
            expected = ν / 2 * (cmax - cmin)
            @test all(parent(total_water) .≈ expected)

            dY = similar(Y)
            exp_tendency! = make_exp_tendency(soil)
            exp_tendency!(dY, Y, p, t0)
            @test dY.soil.∫Fwdt == p.soil.∫Swdz
            @test all(
                parent(p.soil.∫Swdz) .≈ parent(
                    ClimaCore.Fields.zeros(domain.space.surface) .+
                    FT(((cmax - cmin) * -1e-5)),
                ),
            )
            imp_tendency! = make_imp_tendency(soil)
            imp_tendency!(dY, Y, p, t0)
            @test dY.soil.∫Fwdt ==
                  ClimaCore.Fields.zeros(domain.space.surface) .- FT(2.0)
        end
    end
    @testset "Energy Hydrology total energy and water, FT = $FT" begin
        sources = (Soil.PhaseChange{FT}(), FakeSource{FT}())
        params = Soil.EnergyHydrologyParameters(
            FT;
            ν,
            ν_ss_om,
            ν_ss_quartz,
            ν_ss_gravel,
            hydrology_cm = hcm,
            K_sat,
            S_s,
            θ_r,
        )
        boundary_fluxes = (;
            top = WaterHeatBC(;
                water = top_water_flux_bc,
                heat = top_heat_flux_bc,
            ),
            bottom = WaterHeatBC(;
                water = bot_water_flux_bc,
                heat = bot_heat_flux_bc,
            ),
        )
        for domain in domains
            soil = Soil.EnergyHydrology{FT}(;
                parameters = params,
                domain = domain,
                boundary_conditions = boundary_fluxes,
                sources = sources,
            )

            Y, p, cds = initialize(soil)
            ϑ_l0 = ν / 2
            θ_i0 = ν / 5
            T = FT(270.0)
            ρc_s = Soil.volumetric_heat_capacity(
                ϑ_l0,
                θ_i0,
                params.ρc_ds,
                params.earth_param_set,
            )
            ρe_int0 =
                Soil.volumetric_internal_energy.(
                    θ_i0,
                    ρc_s,
                    T,
                    params.earth_param_set,
                )
            Y.soil.ϑ_l .= ϑ_l0
            Y.soil.θ_i .= θ_i0

            Y.soil.ρe_int = ρe_int0
            set_initial_cache! = make_set_initial_cache(soil)
            set_initial_cache!(p, Y, t0)
            total_water = ClimaCore.Fields.zeros(soil.domain.space.surface)
            total_energy = ClimaCore.Fields.zeros(soil.domain.space.surface)
            ClimaLand.total_liq_water_vol_per_area!(total_water, soil, Y, p, t0)
            ρ_ice = LP.ρ_cloud_ice(params.earth_param_set)
            ρ_liq = LP.ρ_cloud_liq(params.earth_param_set)
            expected = (ϑ_l0 + θ_i0 * ρ_ice / ρ_liq) * (cmax - cmin)
            @test all(parent(total_water) .≈ expected)
            ClimaLand.total_energy_per_area!(total_energy, soil, Y, p, t0)
            @test all(parent(total_energy) .≈ ρe_int0 .* (cmax - cmin))

            dY = similar(Y)
            exp_tendency! = make_exp_tendency(soil)
            exp_tendency!(dY, Y, p, t0)
            @test dY.soil.∫Fwdt == p.soil.∫Swdz
            @test dY.soil.∫Fedt == p.soil.∫Sedz
            @test all(
                parent(p.soil.∫Swdz) .≈ parent(
                    ClimaCore.Fields.zeros(domain.space.surface) .+
                    FT(((cmax - cmin) * -1e-5)),
                ),
            )
            imp_tendency! = make_imp_tendency(soil)
            imp_tendency!(dY, Y, p, t0)
            @test dY.soil.∫Fwdt ==
                  ClimaCore.Fields.zeros(domain.space.surface) .- FT(2.0)
            @test dY.soil.∫Fedt ==
                  ClimaCore.Fields.zeros(domain.space.surface) .- FT(2.0)
        end
    end
end

import SciMLBase
using Statistics
import ClimaTimeSteppers as CTS
using ClimaCore
import ClimaParams as CP
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
import ClimaLand
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand.Parameters as LP

FT = Float32;
toml_dict = LP.create_toml_dict(FT)
earth_param_set = LP.LandParameters(toml_dict)

vg_n = FT(2.9)
vg_α = FT(6)
K_sat = FT(4.42 / 3600 / 100) # m/s
hcm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
ν = FT(0.43)
θ_r = FT(0.045)
S_s = FT(1e-3)
ν_ss_om = FT(0.0)
ν_ss_quartz = FT(1.0)
ν_ss_gravel = FT(0.0)

params = ClimaLand.Soil.EnergyHydrologyParameters(
    toml_dict;
    ν,
    ν_ss_om,
    ν_ss_quartz,
    ν_ss_gravel,
    hydrology_cm = hcm,
    K_sat,
    S_s,
    θ_r,
);

surface_water_flux = WaterFluxBC((p, t) -> -K_sat / 10)
bottom_water_flux = WaterFluxBC((p, t) -> 0.0);

surface_heat_flux = HeatFluxBC((p, t) -> 0.0)
bottom_heat_flux = HeatFluxBC((p, t) -> 0.0);

boundary_fluxes = (;
    top = WaterHeatBC(; water = surface_water_flux, heat = surface_heat_flux),
    bottom = WaterHeatBC(; water = bottom_water_flux, heat = bottom_heat_flux),
);

# Domain - single column
zmax = FT(0)
zmin = FT(-0.5)
nelems = 25
soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
z = ClimaCore.Fields.coordinate_field(soil_domain.space.subsurface).z;

# Soil model, and create the prognostic vector Y and cache p:
soil = Soil.EnergyHydrology{FT}(;
    parameters = params,
    domain = soil_domain,
    boundary_conditions = boundary_fluxes,
    sources = (),
)


function hydrostatic_equilibrium(z, z_interface, params)
    (; ν, S_s, hydrology_cm) = params
    (; α, n, m) = hydrology_cm
    if z < z_interface
        return -S_s * (z - z_interface) + ν
    else
        return ν * (1 + (α * (z - z_interface))^n)^(-m)
    end
end
function init_soil!(Y, p, t0, model)
    params = model.parameters
    z = model.domain.fields.z
    FT = eltype(Y.soil.ϑ_l)
    Y.soil.ϑ_l .= hydrostatic_equilibrium.(z, FT(-0.45), params)
    Y.soil.θ_i .= 0
    T = FT(296.15)
    ρc_s = @. Soil.volumetric_heat_capacity(
        Y.soil.ϑ_l,
        FT(0),
        params.ρc_ds,
        params.earth_param_set,
    )
    Y.soil.ρe_int =
        Soil.volumetric_internal_energy.(FT(0), ρc_s, T, params.earth_param_set)
end

# Timestepping:
t0 = Float64(0)
tf = Float64(24 * 3600)
timestepper = CTS.ARS111();
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = 6,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
    ),
);

simulation_dt_small =
    LandSimulation(t0, tf, FT(100), soil; set_ic! = init_soil!)

simulation_dt_large = LandSimulation(
    t0,
    tf,
    FT(900),
    soil;
    set_ic! = init_soil!,
    timestepper = ode_algo,
)

# Solve at small dt
sol_dt_small = solve!(simulation_dt_small);
sol_dt_large = solve!(simulation_dt_large);
@assert sum(isnan.(simulation_dt_large._integrator.u)) == 0

norm(x) = sqrt(mean(parent(x .^ 2)))
rmse(x, y) = sqrt(mean(parent((x .- y) .^ 2)))
@assert rmse(
    simulation_dt_small._integrator.u.soil.ϑ_l,
    simulation_dt_large._integrator.u.soil.ϑ_l,
) / norm(simulation_dt_small._integrator.u.soil.ϑ_l) < sqrt(eps(FT))
@assert rmse(
    simulation_dt_small._integrator.u.soil.ρe_int,
    simulation_dt_large._integrator.u.soil.ρe_int,
) / norm(simulation_dt_small._integrator.u.soil.ρe_int) < sqrt(eps(FT))

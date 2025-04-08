using ClimaComms
using Dates
using Test
import ClimaParams as CP
using ClimaLand
import ClimaLand.Parameters as LP
using ClimaCore
using SciMLBase
using ClimaTimeSteppers
include("full_land_setup.jl")
context = ClimaComms.context()
nelements = (101, 15)
start_date = DateTime(2008)
Δt = 450.0
t0 = 0.0

FT = Float64
earth_param_set = LP.LandParameters(FT)

f_over = FT(3.28) # 1/m
R_sb = FT(1.484e-4 / 1000) # m/s
scalar_soil_params = (; f_over, R_sb)

α_snow = FT(0.67)
scalar_snow_params = (; α_snow, Δt)

# Energy Balance model
ac_canopy = FT(2.5e3)
# Plant Hydraulics and general plant parameters
K_sat_plant = FT(1e-7) # m/s 
ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value is -4 MPa
Weibull_param = FT(4) # unitless, Holtzman's original c param value
a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
plant_ν = FT(1.44e-4)
plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
h_leaf = FT(1.0)
scalar_canopy_params = (;
    ac_canopy,
    K_sat_plant,
    a,
    ψ63,
    Weibull_param,
    plant_ν,
    plant_S_s,
    h_leaf,
)

domain = ClimaLand.global_domain(FT; nelements = nelements)
surface_space = domain.space.surface
start_date = DateTime(2008)
# Forcing data
era5_ncdata_path = joinpath(
    ClimaLand.Artifacts.era5_land_forcing_data2008_folder_path(;
        context,
        lowres = true,
    ),
    "era5_2008_1.0x1.0_lowres.nc",
)
forcing = ClimaLand.prescribed_forcing_era5(
    joinpath(era5_ncdata_path),
    surface_space,
    start_date,
    earth_param_set,
    FT,
)
LAI = ClimaLand.prescribed_lai_modis(
    joinpath(
        ClimaLand.Artifacts.modis_lai_forcing_data_path(; context),
        "Yuan_et_al_2008_1x1.nc",
    ),
    domain.space.surface,
    start_date,
)

land = global_land_model(
    FT,
    scalar_soil_params,
    scalar_canopy_params,
    scalar_snow_params,
    earth_param_set;
    context = context,
    domain = domain,
    forcing = forcing,
    LAI = LAI,
)

Y, p, cds = initialize(land)

# Soil IC
ϑ_l0 = land.soil.parameters.ν ./ 2
θ_i0 = land.soil.parameters.ν ./ 5
T = FT(270.0)
ρc_s = @. Soil.volumetric_heat_capacity(
    ϑ_l0,
    θ_i0,
    land.soil.parameters.ρc_ds,
    earth_param_set,
)
ρe_int0 = @. Soil.volumetric_internal_energy(θ_i0, ρc_s, T, earth_param_set)
Y.soil.ϑ_l .= ϑ_l0
Y.soil.θ_i .= θ_i0

Y.soil.ρe_int = ρe_int0

# Canopy IC
ϑ0 = land.canopy.hydraulics.parameters.ν / 2
CTemp0 = FT(290.5)

Y.canopy.hydraulics.ϑ_l.:1 .= ϑ0

Y.canopy.energy.T .= CTemp0

# Snow IC
S0 = FT(0.5)
S_l0 = FT(0.3)
STemp0 = FT(270)
U0 = Snow.energy_from_T_and_swe(S0, STemp0, land.snow.parameters)

Y.snow.S .= S0
Y.snow.S_l .= S_l0
Y.snow.U .= U0

t0 = 0.0
set_initial_cache! = make_set_initial_cache(land)
set_initial_cache!(p, Y, t0)

# Check total
area_index = p.canopy.hydraulics.area_index.leaf
h_canopy = land.canopy.hydraulics.compartment_surfaces[end]
ρ_ice = LP.ρ_cloud_ice(earth_param_set)
ρ_liq = LP.ρ_cloud_liq(earth_param_set)
int_cache = ClimaCore.Fields.zeros(domain.space.surface)
ClimaCore.Operators.column_integral_definite!(
    int_cache,
    @. (ϑ_l0 + θ_i0 * ρ_ice / ρ_liq)
)
soil_exp = int_cache
canopy_exp = @. area_index * h_canopy .* ϑ0
snow_exp = S0
@test all(parent(p.total_water) .≈ parent(snow_exp .+ canopy_exp .+ soil_exp))

int_cache .*= 0
ClimaCore.Operators.column_integral_definite!(int_cache, ρe_int0)
soil_exp = int_cache
canopy_exp = @. area_index * land.canopy.energy.parameters.ac_canopy * CTemp0
snow_exp = U0
@test all(parent(p.total_energy) .≈ parent(snow_exp .+ canopy_exp .+ soil_exp))


# Check that tendencies conserve
cache = ClimaCore.Fields.zeros(domain.space.surface)
imp_tendency! = ClimaLand.make_imp_tendency(land);
exp_tendency! = ClimaLand.make_exp_tendency(land);
dY_exp = deepcopy(Y) .* 0;
dY_imp = deepcopy(Y) .* 0;
exp_tendency!(dY_exp, Y, p, t0);
p_exp = deepcopy(p)
imp_tendency!(dY_imp, Y, p, t0);

# Canopy
dY_canopy_energy_exp =
    land.canopy.energy.parameters.ac_canopy .*
    p_exp.canopy.hydraulics.area_index.leaf .* dY_exp.canopy.energy.T;
dY_canopy_energy_imp =
    land.canopy.energy.parameters.ac_canopy .*
    p.canopy.hydraulics.area_index.leaf .* dY_imp.canopy.energy.T;
dY_canopy_water_exp =
    p_exp.canopy.hydraulics.area_index.leaf .* dY_exp.canopy.hydraulics.ϑ_l.:1;
dY_canopy_water_imp =
    p.canopy.hydraulics.area_index.leaf .* dY_imp.canopy.hydraulics.ϑ_l.:1;

# Soil
energy_cache_imp = deepcopy(cache)
energy_cache_imp .*= -0
ClimaCore.Operators.column_integral_definite!(
    energy_cache_imp,
    dY_imp.soil.ρe_int,
)
energy_cache_exp = deepcopy(cache)
energy_cache_exp .*= 0
ClimaCore.Operators.column_integral_definite!(
    energy_cache_exp,
    dY_exp.soil.ρe_int,
)

water_cache_imp = deepcopy(cache)
water_cache_imp .*= -0
ClimaCore.Operators.column_integral_definite!(
    water_cache_imp,
    dY_imp.soil.ϑ_l .+ dY_imp.soil.θ_i .* ρ_ice ./ ρ_liq,
)
water_cache_exp = deepcopy(cache)
water_cache_exp .*= 0
ClimaCore.Operators.column_integral_definite!(
    water_cache_exp,
    dY_exp.soil.ϑ_l .+ dY_exp.soil.θ_i .* ρ_ice ./ ρ_liq,
)

# Sum including snow
dY_energy =
    dY_canopy_energy_imp .+ dY_canopy_energy_exp .+ energy_cache_imp .+
    energy_cache_exp .+ dY_imp.snow.U .+ dY_exp.snow.U;
# Expected change based on fluxes at boundaries and sources/sinks in soil column
energy_expected = dY_exp.∫Fedt .+ dY_imp.∫Fedt
@test maximum(
    abs.(
        extrema((energy_expected .- dY_energy) ./ (energy_expected .+ eps(FT)))
    ),
) < sqrt(eps(FT))


dY_water =
    dY_canopy_water_imp .+ dY_canopy_water_exp .+ water_cache_imp .+
    water_cache_exp .+ dY_imp.snow.S .+ dY_exp.snow.S;
# Expected change based on fluxes at boundaries and sources/sinks in soil column
water_expected = dY_exp.∫Fwdt .+ dY_imp.∫Fwdt
@test maximum(
    abs.(extrema((water_expected .- dY_water) ./ (water_expected .+ eps(FT)))),
) < sqrt(eps(FT))

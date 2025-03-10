using ClimaComms
using Dates
using Test
import ClimaParams as CP
using ClimaLand
import ClimaLand.Parameters as LP

context = ClimaComms.context()
nelements = (101, 15)
start_date = DateTime(2008)
Δt = 450.0
t0 = 0.0

FT = Float64
earth_param_set = LP.LandParameters(FT)

f_over = FT(3.28) # 1/m
R_sb = FT(1.484e-4 / 1000) # m/s
other_soil_params = (; f_over, R_sb)

α_snow = FT(0.67)
other_snow_params = (; α_snow,)

# Energy Balance model
ac_canopy = FT(2.5e3)
# Plant Hydraulics and general plant parameters
SAI = FT(0.0) # m2/m2
f_root_to_shoot = FT(3.5)
RAI = FT(1.0)
K_sat_plant = FT(5e-9) # m/s # seems much too small?
ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value is -4 MPa
Weibull_param = FT(4) # unitless, Holtzman's original c param value
a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
plant_ν = FT(1.44e-4)
plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
n_stem = 0
n_leaf = 1
h_stem = FT(0.0)
h_leaf = FT(1.0)
zmax = FT(0.0)
other_canopy_params = (;
    ac_canopy,
    SAI,
    f_root_to_shoot,
    RAI,
    K_sat_plant,
    a,
    ψ63,
    Weibull_param,
    plant_ν,
    plant_S_s,
    n_stem,
    n_leaf,
    h_stem,
    h_leaf,
    zmax,
)

land, Y, p, cds = ClimaLand.land_model_setup(
    FT;
    earth_param_set,
    context,
    nelements,
    start_date,
    t0,
    Δt,
    other_soil_params,
    other_canopy_params,
    other_snow_params,
    era5_lowres = true,
)


# More complex function
set_initial_cache! = make_set_initial_cache(land)
set_initial_cache!(p, Y, t0)
dY = similar(Y)
@. dY = 0

imp_tendency! = make_imp_tendency(land)
imp_tendency!(dY, Y, p, t0)
surface_space = axes(Y.snow.U)
binary_mask = .~parent(surface_space.grid.mask.is_active)[:]
# Test that the masked parts of dY did not update and are still zero
@test extrema(parent(dY.soil.ϑ_l)[:, 1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(dY.soil.θ_i)[:, 1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(dY.soil.ρe_int)[:, 1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(dY.snow.U)[1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(dY.snow.S)[1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(dY.snow.S_l)[1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(dY.canopy.energy.T)[1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(dY.canopy.hydraulics.ϑ_l.:1)[1, 1, 1, binary_mask]) ==
      (0.0, 0.0)

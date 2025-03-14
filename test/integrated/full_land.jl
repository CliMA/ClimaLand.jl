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
scalar_soil_params = (; f_over, R_sb)

α_snow = FT(0.67)
scalar_snow_params = (; α_snow, Δt)

# Energy Balance model
ac_canopy = FT(2.5e3)
# Plant Hydraulics and general plant parameters
K_sat_plant = FT(5e-9) # m/s # seems much too small?
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

domain = ClimaLand.global_domain(FT; nelements = nelements, context = context)
surface_space = domain.space.surface
start_date = DateTime(2008)
# Forcing data
era5_ncdata_path = ClimaLand.Artifacts.era5_land_forcing_data2008_folder_path(;
    context,
    lowres = true,
)
forcing = ClimaLand.prescribed_forcing_era5(
    era5_ncdata_path,
    surface_space,
    start_date,
    earth_param_set,
    FT
)
LAI = ClimaLand.prescribed_lai_modis(
    ClimaLand.Artifacts.modis_lai_forcing_data_path(; context),
    domain.space.surface,
    start_date,
)

land = ClimaLand.global_land_model(
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

@test typeof(land) <: ClimaLand.LandModel

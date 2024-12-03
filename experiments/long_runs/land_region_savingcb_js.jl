# # Global run of land model

# The code sets up and runs the soil/canopy model for 6 hours on a spherical domain,
# using ERA5 data. In this simulation, we have
# turned lateral flow off because horizontal boundary conditions and the
# land/sea mask are not yet supported by ClimaCore.

# Simulation Setup
# Number of spatial elements: 101 in horizontal, 15 in vertical
# Soil depth: 50 m
# Simulation duration: 365 d
# Timestep: 900 s
# Timestepper: ARS343
# Fixed number of iterations: 1
# Jacobian update: every new timestep
# Atmos forcing update: every 3 hours
import SciMLBase
import ClimaComms
ClimaComms.@import_required_backends
import ClimaTimeSteppers as CTS
using ClimaCore
using ClimaUtilities.ClimaArtifacts
import Interpolations
using Insolation

import ClimaDiagnostics
import ClimaAnalysis
import ClimaAnalysis.Visualize as viz
import ClimaUtilities


import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP

using ClimaLand
using ClimaLand.Soil
using ClimaLand.Snow
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP

using Statistics
using CairoMakie
import GeoMakie
using Dates
import NCDatasets

const FT = Float64;

time_interpolation_method = LinearInterpolation(PeriodicCalendar())
regridder_type = :InterpolationsRegridder
context = ClimaComms.context()
device = ClimaComms.device()
t0 = 0.0
tf = 60 * 60.0 * 24 * 450
Δt = 450.0
nelements = (10, 10, 15)

earth_param_set = LP.LandParameters(FT)
depth = FT(50)
center_long, center_lat = FT(-130), FT(57)
delta_m = FT(200_000) #~ few x few degrees
domain = ClimaLand.Domains.HybridBox(;
    xlim = (delta_m, delta_m),
    ylim = (delta_m, delta_m),
    zlim = (-depth, FT(0)),
    nelements = nelements,
    npolynomial = 1,
    longlat = (center_long, center_lat),
    dz_tuple = FT.((10.0, 0.05)),
)
surface_space = domain.space.surface
subsurface_space = domain.space.subsurface

start_date = DateTime(2008)
# Forcing data
era5_artifact_path =
    ClimaLand.Artifacts.era5_land_forcing_data2008_folder_path(; context)
era5_ncdata_path = joinpath(era5_artifact_path, "era5_2008_1.0x1.0.nc")
atmos, radiation = ClimaLand.prescribed_forcing_era5(
    era5_ncdata_path,
    surface_space,
    start_date,
    earth_param_set,
    FT,
)
spatially_varying_soil_params =
    ClimaLand.default_spatially_varying_soil_parameters(
        subsurface_space,
        surface_space,
        FT,
    )
(;
    ν,
    ν_ss_om,
    ν_ss_quartz,
    ν_ss_gravel,
    hydrology_cm,
    K_sat,
    S_s,
    θ_r,
    PAR_albedo_dry,
    NIR_albedo_dry,
    PAR_albedo_wet,
    NIR_albedo_wet,
    f_max,
) = spatially_varying_soil_params
soil_params = Soil.EnergyHydrologyParameters(
    FT;
    ν,
    ν_ss_om,
    ν_ss_quartz,
    ν_ss_gravel,
    hydrology_cm,
    K_sat,
    S_s,
    θ_r,
    PAR_albedo_dry = PAR_albedo_dry,
    NIR_albedo_dry = NIR_albedo_dry,
    PAR_albedo_wet = PAR_albedo_wet,
    NIR_albedo_wet = NIR_albedo_wet,
)
f_over = FT(3.28) # 1/m
R_sb = FT(1.484e-4 / 1000) # m/s
runoff_model = ClimaLand.Soil.Runoff.TOPMODELRunoff{FT}(;
    f_over = f_over,
    f_max = f_max,
    R_sb = R_sb,
)

# Spatially varying canopy parameters from CLM
clm_parameters = ClimaLand.clm_canopy_parameters(surface_space)
(;
    Ω,
    rooting_depth,
    is_c3,
    Vcmax25,
    g1,
    G_Function,
    α_PAR_leaf,
    τ_PAR_leaf,
    α_NIR_leaf,
    τ_NIR_leaf,
) = clm_parameters
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
conductivity_model =
    Canopy.PlantHydraulics.Weibull{FT}(K_sat_plant, ψ63, Weibull_param)
retention_model = Canopy.PlantHydraulics.LinearRetentionCurve{FT}(a)
plant_ν = FT(1.44e-4)
plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
n_stem = 0
n_leaf = 1
h_stem = FT(0.0)
h_leaf = FT(1.0)
zmax = FT(0.0)
h_canopy = h_stem + h_leaf
compartment_midpoints =
    n_stem > 0 ? [h_stem / 2, h_stem + h_leaf / 2] : [h_leaf / 2]
compartment_surfaces =
    n_stem > 0 ? [zmax, h_stem, h_canopy] : [zmax, h_leaf]

z0_m = FT(0.13) * h_canopy
z0_b = FT(0.1) * z0_m


soilco2_ps = Soil.Biogeochemistry.SoilCO2ModelParameters(FT)

soil_args = (domain = domain, parameters = soil_params)
soil_model_type = Soil.EnergyHydrology{FT}

# Soil microbes model
soilco2_type = Soil.Biogeochemistry.SoilCO2Model{FT}

# soil microbes args
Csom = ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))

# Set the soil CO2 BC to being atmospheric CO2
soilco2_top_bc = Soil.Biogeochemistry.AtmosCO2StateBC()
soilco2_bot_bc = Soil.Biogeochemistry.SoilCO2FluxBC((p, t) -> 0.0) # no flux
soilco2_sources = (Soil.Biogeochemistry.MicrobeProduction{FT}(),)

soilco2_boundary_conditions = (; top = soilco2_top_bc, bottom = soilco2_bot_bc)

soilco2_args = (;
    boundary_conditions = soilco2_boundary_conditions,
    sources = soilco2_sources,
    domain = domain,
    parameters = soilco2_ps,
)

# Now we set up the canopy model, which we set up by component:
# Component Types
canopy_component_types = (;
    autotrophic_respiration = Canopy.AutotrophicRespirationModel{FT},
    radiative_transfer = Canopy.TwoStreamModel{FT},
    photosynthesis = Canopy.FarquharModel{FT},
    conductance = Canopy.MedlynConductanceModel{FT},
    hydraulics = Canopy.PlantHydraulicsModel{FT},
    energy = Canopy.BigLeafEnergyModel{FT},
)
# Individual Component arguments
# Set up autotrophic respiration
autotrophic_respiration_args =
    (; parameters = Canopy.AutotrophicRespirationParameters(FT))
# Set up radiative transfer
radiative_transfer_args = (;
    parameters = Canopy.TwoStreamParameters(
        FT;
        Ω,
        α_PAR_leaf,
        τ_PAR_leaf,
        α_NIR_leaf,
        τ_NIR_leaf,
        G_Function,
    )
)
# Set up conductance
conductance_args = (; parameters = Canopy.MedlynConductanceParameters(FT; g1))
# Set up photosynthesis
photosynthesis_args =
    (; parameters = Canopy.FarquharParameters(FT, is_c3; Vcmax25 = Vcmax25))
# Set up plant hydraulics
era5_lai_artifact_path = ClimaLand.Artifacts.era5_lai_forcing_data2008_folder_path(; context)
era5_lai_ncdata_path = joinpath(era5_lai_artifact_path, "era5_2008_1.0x1.0_lai.nc")
LAIfunction = ClimaLand.prescribed_lai_era5(
    era5_lai_ncdata_path,
    surface_space,
    start_date;
    time_interpolation_method = time_interpolation_method,
)
ai_parameterization = Canopy.PrescribedSiteAreaIndex{FT}(LAIfunction, SAI, RAI)

plant_hydraulics_ps = Canopy.PlantHydraulics.PlantHydraulicsParameters(;
    ai_parameterization = ai_parameterization,
    ν = plant_ν,
    S_s = plant_S_s,
    rooting_depth = rooting_depth,
    conductivity_model = conductivity_model,
    retention_model = retention_model,
)
plant_hydraulics_args = (
    parameters = plant_hydraulics_ps,
    n_stem = n_stem,
    n_leaf = n_leaf,
    compartment_midpoints = compartment_midpoints,
    compartment_surfaces = compartment_surfaces,
)

energy_args = (parameters = Canopy.BigLeafEnergyParameters{FT}(ac_canopy),)

# Canopy component args
canopy_component_args = (;
    autotrophic_respiration = autotrophic_respiration_args,
    radiative_transfer = radiative_transfer_args,
    photosynthesis = photosynthesis_args,
    conductance = conductance_args,
    hydraulics = plant_hydraulics_args,
    energy = energy_args,
)

# Other info needed
shared_params = Canopy.SharedCanopyParameters{FT, typeof(earth_param_set)}(
    z0_m,
    z0_b,
    earth_param_set,
)

canopy_model_args = (;
    parameters = shared_params,
    domain = ClimaLand.obtain_surface_domain(domain),
)
# Snow model
snow_parameters = SnowParameters{FT}(
    Δt;
    earth_param_set = earth_param_set,
)
snow_args = (;
    parameters = snow_parameters,
    domain = ClimaLand.obtain_surface_domain(domain),
)
snow_model_type = Snow.SnowModel
# Integrated plant hydraulics and soil model
land_input = (
    atmos = atmos,
    radiation = radiation,
    runoff = runoff_model,
    soil_organic_carbon = Csom,
)
land = LandModel{FT}(;
    soilco2_type = soilco2_type,
    soilco2_args = soilco2_args,
    land_args = land_input,
    soil_model_type = soil_model_type,
    soil_args = soil_args,
    canopy_component_types = canopy_component_types,
    canopy_component_args = canopy_component_args,
    canopy_model_args = canopy_model_args,
    snow_args = snow_args,
    snow_model_type = snow_model_type,
)

Y, p, cds = initialize(land)
Y.snow.S .= 0.0
Y.snow.U .= 0.0
init_soil(ν, θ_r) = θ_r + (ν - θ_r) / 2
Y.soil.ϑ_l .= init_soil.(ν, θ_r)
Y.soil.θ_i .= FT(0.0)
T = FT(276.85)
ρc_s =
    Soil.volumetric_heat_capacity.(
        Y.soil.ϑ_l,
        Y.soil.θ_i,
        soil_params.ρc_ds,
        soil_params.earth_param_set,
    )
Y.soil.ρe_int .=
    Soil.volumetric_internal_energy.(
        Y.soil.θ_i,
        ρc_s,
        T,
        soil_params.earth_param_set,
    )
Y.soilco2.C .= FT(0.000412) # set to atmospheric co2, mol co2 per mol air
Y.canopy.hydraulics.ϑ_l.:1 .= plant_ν
evaluate!(Y.canopy.energy.T, atmos.T, t0)

set_initial_cache! = make_set_initial_cache(land)
exp_tendency! = make_exp_tendency(land)
imp_tendency! = ClimaLand.make_imp_tendency(land)
jacobian! = ClimaLand.make_jacobian(land)
set_initial_cache!(p, Y, t0);

# set up jacobian info
jac_kwargs = (; jac_prototype = ImplicitEquationJacobian(Y), Wfact = jacobian!)

prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction(
        T_exp! = exp_tendency!,
        T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
        dss! = ClimaLand.dss!,
    ),
    Y,
    (t0, tf),
    p,
)

updateat = Array(t0:(3600 * 3):tf)
drivers = ClimaLand.get_drivers(land)
updatefunc = ClimaLand.make_update_drivers(drivers)
t_first_save = t0 + 395*24*3600.0
dt_save = Δt
saveat = Array(t_first_save:dt_save:tf)
sv = (;
    t = Array{Float64}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
)
saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat)

driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
cb = SciMLBase.CallbackSet(driver_cb, saving_cb)


# Define timestepper and ODE algorithm
stepper = CTS.ARS111()
ode_algo = CTS.IMEXAlgorithm(
    stepper,
    CTS.NewtonsMethod(
        max_iters = 3,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
    ),
)
sol = SciMLBase.solve(
    prob,
    ode_algo;
    dt = Δt,
    callback = cb,
    adaptive = false,
    saveat = saveat,
);
extract2d(x) = mean(Array(parent(x))[1,1,1,26])
extract3d(x) =
    mean(Array(parent(ClimaLand.top_center_to_surface(x)))[1, 1, 1, 26])
id_last = length(sol.t)
S = [
    extract3d((sv.saveval[k].soil.θ_l .- θ_r) ./ (ν .- θ_r)) for k in 1:id_last
];
ϑ_l = [extract3d(sol.u[k].soil.ϑ_l) for k in 1:id_last];
SWE = [extract2d(sol.u[k].snow.S) for k in 1:id_last];
T_snow = [extract2d(sv.saveval[k].snow.T) for k in 1:id_last];
Δz_soil  = [extract2d(sv.saveval[k].effective_soil_sfc_depth) for k in 1:id_last];
T_eff_soil  = [extract2d(sv.saveval[k].effective_soil_sfc_T) for k in 1:id_last];

Δz_snow  = [extract2d(sv.saveval[k].snow.z) for k in 1:id_last];
ghf = [extract2d(sv.saveval[k].ground_heat_flux) for k in 1:id_last];
    κ_snow = [extract2d(sv.saveval[k].snow.κ  ) for k in 1:id_last];
    κ_soil = [extract3d(sv.saveval[k].soil.κ  ) for k in 1:id_last];
ψ = [extract3d(sv.saveval[k].soil.ψ) for k in 1:id_last];
T_soil = [extract3d(sv.saveval[k].soil.T) for k in 1:id_last];

θ_i = [extract3d(sol.u[k].soil.θ_i) for k in 1:id_last];
top_bc = [extract2d(sv.saveval[k].soil.top_bc.water) for k in 1:id_last]
top_bc_heat = [extract2d(sv.saveval[k].soil.top_bc.heat) for k in 1:id_last]

precip = [extract2d(sv.saveval[k].drivers.P_liq) for k in 1:id_last]
evap = [
    extract2d(sv.saveval[k].soil.turbulent_fluxes.vapor_flux_liq) for
    k in 1:id_last
]
sub = [
    extract2d(sv.saveval[k].soil.turbulent_fluxes.vapor_flux_ice) for
    k in 1:id_last
]
extraction = [extract3d(sv.saveval[k].root_extraction) for k in 1:id_last]
T_canopy = [extract2d(sol.u[k].canopy.energy.T) for k in 1:id_last]

T_air = [extract2d(sv.saveval[k].drivers.T) for k in 1:id_last]

days = sv.t[1:id_last] ./ 24 ./ 3600


output_dir = joinpath(pwd(), "land_region_savingcb_output_11-26")
!isdir(output_dir) && mkdir(output_dir)

using DelimitedFiles

# Get all state and cache variables we accessed above
vardict = ("S" => S,
    "ϑ_l" => ϑ_l,
    "SWE" => SWE,
    "T_snow" => T_snow,
    "Δz_soil" => Δz_soil,
    "T_eff_soil" => T_eff_soil,
    "Δz_snow" => Δz_snow,
    "ghf" => ghf,
    "ψ" => ψ,
    "T_soil" => T_soil,
    "θ_i" => θ_i,
    "top_bc" => top_bc,
    "top_bc_heat" => top_bc_heat,
    "precip" => precip,
    "evap" => evap,
    "sub" => sub,
    "extraction" => extraction,
    "T_canopy" => T_canopy,
    "T_air" => T_air)

# Save everything to comma-delimited text files
for (varname, var) in vardict
    open(
        joinpath(output_dir, "$varname.txt"),
        "w",
    ) do io
        writedlm(io, var, ',')
    end;
end

fig = CairoMakie.Figure(size = (1500, 1500))
ax = CairoMakie.Axis(fig[1, 1], xlabel = "Time", ylabel = "Water Content")
#lines!(ax, days, ϑ_l .- extract3d(θ_r), label = "ϑ_l- θ_r")
lines!(ax, days, θ_i, label = "θ_i")
lines!(ax, days, S, label = "S_l")
axislegend(ax)



ax = CairoMakie.Axis(
    fig[2, 2],
    xlabel = "Time",
    ylabel = "Matric Potential",
    yscale = log10,
)
lines!(ax, days, abs.(ψ), label = "ψ(m)")
axislegend(ax)
ax = CairoMakie.Axis(fig[2, 3], xlabel = "Time", ylabel = "Temperature")
lines!(ax, days, T, label = "T soil")
lines!(ax, days, T_canopy, label = "T canopy")
lines!(ax, days, T_eff_soil, label = "T eff_soil")
lines!(ax, days, T_air, label = "T air")
lines!(ax, days, T_snow, label = "T snow")
axislegend(ax)
ax = CairoMakie.Axis(fig[2, 1], xlabel = "Time", ylabel = "Water Fluxes")
lines!(ax, days, evap, label = "E")
lines!(ax, days, sub, label = "S")
lines!(ax, days, precip, label = "precip")
lines!(ax, days, top_bc, label = "BC")
lines!(ax, days, extraction ./ 0.025, label = "extraction")
axislegend(ax)

ax = CairoMakie.Axis(fig[1, 2], xlabel = "Time", ylabel = "Heat Fluxes")
lines!(ax, days, ghf, label = "ghf")
lines!(ax, days, top_bc_heat, label = "BC")
axislegend(ax)
ax = CairoMakie.Axis(fig[1, 3], xlabel = "Time", ylabel = "SWE")
lines!(ax, days, SWE, label = "SWE")
axislegend(ax)

canopy_water = [extract2d(sol.u[k].canopy.hydraulics.ϑ_l .- plant_ν) for k in 1:id_last]
ax = CairoMakie.Axis(fig[1, 4], xlabel = "Time", ylabel = "Canopy Moisture")
lines!(ax, days, canopy_water, label = "canopy_water")
axislegend(ax)
Tra = [extract2d(sv.saveval[k].canopy.energy.turbulent_fluxes.transpiration) for k in 1:id_last]

ax = CairoMakie.Axis(fig[2,4], xlabel = "Time", ylabel = "Canopy water fluxes")
lines!(ax, days, Tra, label = "Tra")
axislegend(ax)
CairoMakie.save("snowy.png", fig)

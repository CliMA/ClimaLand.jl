"""
This experiment tests running the Ozark site (US-MOz) using plant parameters
defined by plant functional types instead of fully site-specific parameters.
"""

import ClimaLand

import SciMLBase
import ClimaTimeSteppers as CTS
using ClimaCore
import ClimaParams as CP
using CairoMakie
using Statistics
using Dates
using Insolation
using StatsBase

using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
import ClimaLand
import ClimaLand.Parameters as LP
import ClimaUtilities.OutputPathGenerator: generate_output_path
using ClimaDiagnostics
using ClimaUtilities
using DelimitedFiles
FluxnetSimulationsExt =
    Base.get_extension(ClimaLand, :FluxnetSimulationsExt).FluxnetSimulationsExt;

const FT = Float64
earth_param_set = LP.LandParameters(FT)
climaland_dir = pkgdir(ClimaLand)
include(joinpath(climaland_dir, "experiments/integrated/fluxnet/plot_utils.jl"))

site_ID = "US-MOz"

# Read all site-specific domain parameters from the simulation file for the site
include(
    joinpath(
        climaland_dir,
        "experiments/integrated/fluxnet/$site_ID/$(site_ID)_simulation.jl",
    ),
)

include(
    joinpath(climaland_dir, "experiments/integrated/fluxnet/fluxnet_domain.jl"),
)

# Read in the site-specific parameters for all parameters not defined by the PFTs
include(
    joinpath(
        climaland_dir,
        "experiments/integrated/fluxnet/$site_ID/$(site_ID)_parameters.jl",
    ),
)

land_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = dz_tuple,
    longlat = (long, lat),
)
canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)
prognostic_land_components = (:canopy, :soil, :soilco2)

# Define the PFT land cover percentages for the Ozark site. Currently we only
# use the dominant PFT, which for Ozark is deciduous broadleaf temperate trees.
pft_pcts = [
    0.0, # NET_Temp
    0.0, # NET_Bor
    0.0, # NDT_Bor
    0.0, # BET_Trop
    0.0, # BET_Temp
    0.0, # BDT_Trop
    1.0, # BDT_Temp
    0.0, # BDT_Bor
    0.0, # BES_Temp
    0.0, # BDS_Temp
    0.0, # BDT_Bor
    0.0, # C3G_A
    0.0, # C3G_NA
    0.0, # C4G
]

# Load the PFT parameters into the namespace
(
    Ω,
    α_PAR_leaf,
    α_NIR_leaf,
    τ_PAR_leaf,
    τ_NIR_leaf,
    ϵ_canopy,
    χl,
    ac_canopy,
    g1,
    Vcmax25,
    f_root_to_shoot,
    K_sat_plant,
    ψ63,
    plant_ν,
    rooting_depth,
) = FT.(params_from_pfts(pft_pcts))

# For now, replace ac_canopy with larger value in order to increase
# stability of timestepping
ac_canopy = ac_canopy * 3

# This sets up the timestepper for the simulation
include(
    joinpath(
        climaland_dir,
        "experiments/integrated/fluxnet/fluxnet_simulation.jl",
    ),
)
start_date = DateTime(2010) + Hour(time_offset)
(; atmos, radiation) = FluxnetSimulationsExt.prescribed_forcing_fluxnet(
    site_ID,
    lat,
    long,
    time_offset,
    atmos_h,
    start_date,
    earth_param_set,
    FT,
)


# Now we set up the model. For the soil model, we pick
# a model type and model args:
soil_domain = land_domain
soil_albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(;
    PAR_albedo = soil_α_PAR,
    NIR_albedo = soil_α_NIR,
)

runoff = ClimaLand.Soil.Runoff.SurfaceRunoff()
retention_parameters = (;
    ν = soil_ν,
    θ_r,
    K_sat = soil_K_sat,
    hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n),
)
composition_parameters = (; ν_ss_om, ν_ss_quartz, ν_ss_gravel)
soil_forcing = (; atmos, radiation)
soil = Soil.EnergyHydrology{FT}(
    soil_domain,
    soil_forcing,
    earth_param_set;
    prognostic_land_components,
    additional_sources = (ClimaLand.RootExtraction{FT}(),),
    albedo = soil_albedo,
    runoff,
    retention_parameters,
    composition_parameters,
    S_s = soil_S_s,
    z_0m = z_0m_soil,
    z_0b = z_0b_soil,
    emissivity = soil_ϵ,
)

# Soil microbes model
soil_organic_carbon =
    ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))
co2_prognostic_soil = Soil.Biogeochemistry.PrognosticMet(soil.parameters)
drivers = Soil.Biogeochemistry.SoilDrivers(
    co2_prognostic_soil,
    soil_organic_carbon,
    atmos,
)
soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(; domain = soil_domain, drivers)

# Canopy model - Set up individual Component arguments with non-default parameters
# Set up radiative transfer
radiation_parameters = (;
    Ω,
    G_Function = CLMGFunction(χl),
    α_PAR_leaf,
    τ_PAR_leaf,
    α_NIR_leaf,
    τ_NIR_leaf,
)
radiative_transfer =
    Canopy.TwoStreamModel{FT}(canopy_domain; radiation_parameters, ϵ_canopy)

# Set up conductance
conductance = Canopy.MedlynConductanceModel{FT}(canopy_domain; g1)

# Set up photosynthesis
photosynthesis_parameters = (; is_c3 = FT(1), Vcmax25)
photosynthesis =
    FarquharModel{FT}(canopy_domain; photosynthesis_parameters, sc, pc)

# Set up plant hydraulics
(; LAI, maxLAI) =
    FluxnetSimulationsExt.prescribed_LAI_fluxnet(site_ID, start_date)
RAI = maxLAI * f_root_to_shoot
hydraulics = Canopy.PlantHydraulicsModel{FT}(
    canopy_domain,
    LAI;
    n_stem,
    n_leaf,
    h_stem,
    h_leaf,
    SAI,
    RAI,
    ν = plant_ν,
    S_s = plant_S_s,
    conductivity_model,
    retention_model,
    rooting_depth,
)

# Set up energy model
energy = Canopy.BigLeafEnergyModel{FT}(; ac_canopy)

ground = ClimaLand.PrognosticSoilConditions{FT}()
canopy_forcing = (; atmos, radiation, ground)
# Construct the canopy model using defaults for autotrophic respiration and SIF models
canopy = Canopy.CanopyModel{FT}(
    canopy_domain,
    canopy_forcing,
    LAI,
    earth_param_set;
    z_0m = z0_m,
    z_0b = z0_b,
    prognostic_land_components,
    radiative_transfer,
    photosynthesis,
    conductance,
    hydraulics,
    energy,
)

# Integrated plant hydraulics and soil model
land = SoilCanopyModel{FT}(soilco2, soil, canopy)

Y, p, cds = initialize(land)
exp_tendency! = make_exp_tendency(land)
imp_tendency! = make_imp_tendency(land);
jacobian! = make_jacobian(land);
jac_kwargs =
    (; jac_prototype = ClimaLand.FieldMatrixWithSolver(Y), Wfact = jacobian!);

FluxnetSimulationsExt.set_fluxnet_ic!(Y, site_ID, start_date, time_offset, land)
set_initial_cache! = make_set_initial_cache(land)
set_initial_cache!(p, Y, t0);

# Callbacks
outdir = joinpath(pkgdir(ClimaLand), "experiments/integrated/fluxnet/out")
output_dir = ClimaUtilities.OutputPathGenerator.generate_output_path(outdir)

output_writer = ClimaDiagnostics.Writers.DictWriter()

short_names_1D = [
    "sif",
    "ra",
    "gs",
    "gpp",
    "ct",
    "swu",
    "lwu",
    "er",
    "et",
    "msf",
    "shf",
    "lhf",
    "rn",
]
short_names_2D = ["swc", "tsoil", "si"]
output_vars = [short_names_1D..., short_names_2D...]

diags = ClimaLand.default_diagnostics(
    land,
    start_date;
    output_writer = output_writer,
    output_vars,
    average_period = :hourly,
)

diagnostic_handler =
    ClimaDiagnostics.DiagnosticsHandler(diags, Y, p, t0, dt = dt);

diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler);

## How often we want to update the drivers. Note that this uses the defined `t0` and `tf`
## defined in the simulatons file
data_dt = Float64(FluxnetSimulationsExt.get_data_dt(site_ID))
updateat = Array(t0:data_dt:tf)
model_drivers = ClimaLand.get_drivers(land)
updatefunc = ClimaLand.make_update_drivers(model_drivers)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
cb = SciMLBase.CallbackSet(driver_cb, diag_cb)


prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction(
        T_exp! = exp_tendency!,
        T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
        dss! = ClimaLand.dss!,
    ),
    Y,
    (t0, tf),
    p,
);

sol = SciMLBase.solve(prob, ode_algo; dt = dt, callback = cb);

ClimaLand.Diagnostics.close_output_writers(diags)

# Extract model output from the saved diagnostics
hourly_diag_name = short_names_1D .* "_1h_average"
hourly_diag_name_2D = short_names_2D .* "_1h_average"


# diagnostic_as_vectors()[2] is a vector of a variable,
# whereas diagnostic_as_vectors()[1] is a vector or time associated with that variable.
# We index to only extract the period post-spinup.
SIF, AR, g_stomata, GPP, canopy_T, SW_u, LW_u, ER, ET, β, SHF, LHF, Rn = [
    ClimaLand.Diagnostics.diagnostic_as_vectors(output_writer, diag_name)[2][(N_spinup_days * 24):end]
    for diag_name in hourly_diag_name
]

swc, soil_T, si = [
    ClimaLand.Diagnostics.diagnostic_as_vectors(
        output_writer,
        diag_name;
        layer = 20, #surface layer
    )[2][(N_spinup_days * 24):end] for diag_name in hourly_diag_name_2D
]
dt_save = 3600.0 # hourly diagnostics
# Number of days to plot post spinup
num_days = N_days - N_spinup_days
model_times = Array(0:dt_save:(num_days * S_PER_DAY)) .+ t_spinup # post spin-up

# convert units for GPP and ET
GPP = GPP .* 1e6 # mol to μmol
AR = AR .* 1e6
ET = ET .* 24 .* 3600

# For the data, we also restrict to post-spinup period
data_id_post_spinup = Array(Int64(t_spinup ÷ data_dt):1:Int64(tf ÷ data_dt))
data_times = Array(0:data_dt:(num_days * S_PER_DAY)) .+ t_spinup
# Plotting
savedir =
    generate_output_path("experiments/integrated/fluxnet/$site_ID/out/pft/")

if !isdir(savedir)
    mkdir(savedir)
end

# Plot model diurnal cycles without data comparisons
comparison_data = FluxnetSimulationsExt.get_comparison_data(
    site_ID,
    time_offset,
    start_date,
    FT,
)
if comparison_data.GPP.absent
    plot_daily_avg(
        "GPP",
        GPP,
        dt_save,
        num_days,
        "μmol/m^2/s",
        savedir,
        "Model",
    )
else
    GPP_data = comparison_data.GPP.values[data_id_post_spinup] .* 1e6
    plot_avg_comp(
        "GPP",
        GPP,
        dt_save,
        GPP_data,
        FT(data_dt),
        num_days,
        "μmol/m^2/s",
        savedir,
    )
end

# Upwelling shortwave radiation is referred to as outgoing in the data
if comparison_data.SW_u.absent
    plot_daily_avg("SW up", SW_u, dt_save, num_days, "w/m^2", savedir, "model")
else
    SW_u_data = comparison_data.SW_u.values[data_id_post_spinup]
    plot_avg_comp(
        "SW up",
        SW_u,
        dt_save,
        SW_u_data,
        FT(data_dt),
        num_days,
        "W/m^2",
        savedir,
    )
end

# Upwelling longwave radiation is referred to outgoing in the data
if comparison_data.LW_u.absent
    plot_daily_avg("LW up", LW_u, dt_save, num_days, "w/m^2", savedir, "model")
else
    LW_u_data = comparison_data.LW_u.values[data_id_post_spinup]
    plot_avg_comp(
        "LW up",
        LW_u,
        dt_save,
        LW_u_data,
        FT(data_dt),
        num_days,
        "W/m^2",
        savedir,
    )
end

# ET
if comparison_data.LE.absent
    plot_daily_avg("ET", ET, dt_save, num_days, "mm/day", savedir, "Model")
else
    measured_T =
        comparison_data.LE.values ./ (LP.LH_v0(earth_param_set) * 1000) .*
        (1e3 * 24 * 3600)
    ET_data = measured_T[data_id_post_spinup]
    plot_avg_comp(
        "ET",
        ET,
        dt_save,
        ET_data,
        FT(data_dt),
        num_days,
        "mm/day",
        savedir,
    )
end

# Sensible Heat Flux
if comparison_data.H.absent
    plot_daily_avg("SHF", SHF, dt_save, num_days, "w/m^2", savedir, "Model")
else
    SHF_data = comparison_data.H.values[data_id_post_spinup]
    plot_avg_comp(
        "SHF",
        SHF,
        dt_save,
        SHF_data,
        FT(data_dt),
        num_days,
        "W/m^2",
        savedir,
    )
end

# Latent Heat Flux
if comparison_data.LE.absent
    plot_daily_avg("LHF", LHF, dt_save, num_days, "w/m^2", savedir)
else
    LHF_data = comparison_data.LE.values[data_id_post_spinup]
    plot_avg_comp(
        "LHF",
        LHF,
        dt_save,
        LHF_data,
        FT(data_dt),
        num_days,
        "W/m^2",
        savedir,
    )
end

# Water stress factor
plot_daily_avg("moisture_stress", β, dt_save, num_days, "", savedir, "Model")

# Stomatal conductance
plot_daily_avg(
    "stomatal_conductance",
    g_stomata,
    dt_save,
    num_days,
    "m s^-1",
    savedir,
    "Model",
)

if isfile(
    joinpath(
        climaland_dir,
        "experiments/integrated/fluxnet/$site_ID/Artifacts.toml",
    ),
)
    rm(
        joinpath(
            climaland_dir,
            "experiments/integrated/fluxnet/$site_ID/Artifacts.toml",
        ),
    )
end

# Water content in soil and snow
# Soil water content
# Current resolution has the first layer at 0.1 cm, the second at 5cm.
fig = Figure(size = (1000, 800), fontsize = 20)
ax1 = Axis(fig[2, 1], xlabel = "Days", ylabel = "SWC [m/m]")
limits!(
    ax1,
    minimum(model_times ./ 3600 ./ 24),
    maximum(model_times ./ 3600 ./ 24),
    0.05,
    0.6,
)
lines!(ax1, model_times ./ 3600 ./ 24, swc, label = "1.25cm", color = "blue")
lines!(
    ax1,
    model_times ./ 3600 ./ 24,
    si,
    color = "cyan",
    label = "Ice, 1.25cm",
)

if !comparison_data.SWC.absent
    lines!(
        ax1,
        data_times ./ 3600 ./ 24,
        comparison_data.SWC.values[data_id_post_spinup],
        label = "Data",
    )
end
ax3 = Axis(
    fig[1, 1],
    xlabel = "",
    ylabel = "Precipitation [mm/day]",
    xticklabelsvisible = false,
)
limits!(
    ax3,
    minimum(model_times ./ 3600 ./ 24),
    maximum(model_times ./ 3600 ./ 24),
    -500,
    0.0,
)

lines!(
    ax3,
    data_times ./ 3600 ./ 24,
    (comparison_data.P.values .* (-1e3 * 24 * 3600))[data_id_post_spinup],
    label = "Total precip (data)",
)

CairoMakie.save(joinpath(savedir, "ground_water_content.png"), fig)



fig2 = Figure(size = (1500, 800))
ax12 = Axis(fig2[1, 1], xlabel = "Days", ylabel = "Temperature (K)")
limits!(
    ax12,
    minimum(model_times ./ 3600 ./ 24),
    maximum(model_times ./ 3600 ./ 24),
    265,
    315,
)

lines!(
    ax12,
    model_times ./ 3600 ./ 24,
    soil_T,
    label = "1.25cm",
    color = "blue",
)

if !comparison_data.TS.absent
    lines!(
        ax12,
        data_times ./ 3600 ./ 24,
        comparison_data.TS.values[data_id_post_spinup],
        label = "Data",
    )
end

CairoMakie.save(joinpath(savedir, "soil_temperature.png"), fig2)

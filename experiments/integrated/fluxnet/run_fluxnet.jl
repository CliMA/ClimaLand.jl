import SciMLBase
import ClimaTimeSteppers as CTS
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
import ClimaParams as CP
using Statistics
using Dates
using Insolation
using StatsBase

using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Snow
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

# Read in the site to be run from the command line
if length(ARGS) < 1
    error("Must provide site ID on command line")
end

site_ID = ARGS[1]

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

# Read all site-specific parameters from the parameter file for the site
include(
    joinpath(
        climaland_dir,
        "experiments/integrated/fluxnet/$site_ID/$(site_ID)_parameters.jl",
    ),
)

# This reads in the data from the flux tower site and creates
# the atmospheric and radiative driver structs for the model
include(
    joinpath(
        climaland_dir,
        "experiments/integrated/fluxnet/fluxnet_simulation.jl",
    ),
)
(start_date, end_date) =
    FluxnetSimulationsExt.get_data_dates(site_ID, time_offset)
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
(; LAI, maxLAI) =
    FluxnetSimulationsExt.prescribed_LAI_fluxnet(site_ID, start_date)
RAI = maxLAI * f_root_to_shoot
capacity = plant_ν * maxLAI * h_leaf * FT(1000)

# Now we set up the model. For the soil model, we pick
# a model type and model args:
soil_domain = land_domain
soil_albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(;
    PAR_albedo = soil_α_PAR,
    NIR_albedo = soil_α_NIR,
)
soil_ps = Soil.EnergyHydrologyParameters(
    FT;
    ν = soil_ν,
    ν_ss_om,
    ν_ss_quartz,
    ν_ss_gravel,
    hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n),
    K_sat = soil_K_sat,
    S_s = soil_S_s,
    θ_r,
    z_0m = z_0m_soil,
    z_0b = z_0b_soil,
    emissivity = soil_ϵ,
    albedo = soil_albedo,
);

soil_args = (domain = soil_domain, parameters = soil_ps)
soil_model_type = Soil.EnergyHydrology{FT}

# Soil microbes model
soilco2_type = Soil.Biogeochemistry.SoilCO2Model{FT}
soilco2_ps = Soil.Biogeochemistry.SoilCO2ModelParameters(FT)
Csom = ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))
soilco2_args = (; domain = soil_domain, parameters = soilco2_ps)

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
    (; parameters = AutotrophicRespirationParameters(FT))
# Set up radiative transfer
radiative_transfer_args = (;
    parameters = TwoStreamParameters(
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
conductance_args = (; parameters = MedlynConductanceParameters(FT; g1))
# Set up photosynthesis
is_c3 = FT(1) # set the photosynthesis mechanism to C3
photosynthesis_args =
    (; parameters = FarquharParameters(FT, is_c3; Vcmax25 = Vcmax25))
# Set up plant hydraulics
ai_parameterization = PrescribedSiteAreaIndex{FT}(LAI, SAI, RAI)

plant_hydraulics_ps = PlantHydraulics.PlantHydraulicsParameters(;
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
shared_params = SharedCanopyParameters{FT, typeof(earth_param_set)}(
    z0_m,
    z0_b,
    earth_param_set,
)

canopy_model_args = (; parameters = shared_params, domain = canopy_domain)

# Snow model
snow_parameters = SnowParameters{FT}(dt; earth_param_set = earth_param_set);
snow_args = (; parameters = snow_parameters, domain = canopy_domain);
snow_model_type = Snow.SnowModel
# Integrated plant hydraulics and soil model
land_input = (
    atmos = atmos,
    radiation = radiation,
    soil_organic_carbon = Csom,
    runoff = ClimaLand.Soil.Runoff.SurfaceRunoff(),
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

FluxnetSimulationsExt.set_fluxnet_ic!(Y, site_ID, start_date, time_offset, land)
set_initial_cache! = make_set_initial_cache(land)
set_initial_cache!(p, Y, t0);

exp_tendency! = make_exp_tendency(land)
imp_tendency! = make_imp_tendency(land)
jacobian! = make_jacobian(land);
jac_kwargs =
    (; jac_prototype = ClimaLand.FieldMatrixWithSolver(Y), Wfact = jacobian!);


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
    "swe",
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

## How often we want to update the drivers. Note that this uses the defined `t0`, and `tf`
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

@time sol = SciMLBase.solve(prob, ode_algo; dt = dt, callback = cb);

ClimaLand.Diagnostics.close_output_writers(diags)
hourly_diag_name = short_names_1D .* "_1h_average"
hourly_diag_name_2D = short_names_2D .* "_1h_average"


# diagnostic_as_vectors()[2] is a vector of a variable,
# whereas diagnostic_as_vectors()[1] is a vector or time associated with that variable.
# We index to only extract the period post-spinup.
SIF, AR, g_stomata, GPP, canopy_T, SW_u, LW_u, ER, ET, β, SHF, LHF, Rn, SWE = [
    ClimaLand.Diagnostics.diagnostic_as_vectors(output_writer, diag_name)[2][(N_spinup_days * 24):end]
    for diag_name in hourly_diag_name
]

swc, soil_T, si = [
    ClimaLand.Diagnostics.diagnostic_as_vectors(
        output_writer,
        diag_name;
        layer = nelements, #surface layer
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
savedir = generate_output_path("experiments/integrated/fluxnet/$site_ID/out/")

if !isdir(savedir)
    mkdir(savedir)
end

# Plot model diurnal cycles without data comparisons
# Autotrophic Respiration
plot_daily_avg(
    "AutoResp",
    AR,
    dt_save,
    num_days,
    "μmol/m^2/s",
    savedir,
    "Model",
)

# Plot all comparisons of model diurnal cycles to data diurnal cycles
# GPP
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
fig = Figure(size = (1500, 800), fontsize = 20)
ax1 = Axis(fig[3, 1], xlabel = "Days", ylabel = "SWC [m/m]")
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
ax2 =
    Axis(fig[2, 1], xlabel = "", ylabel = "SWE [m]", xticklabelsvisible = false)
limits!(
    ax2,
    minimum(model_times ./ 3600 ./ 24),
    maximum(model_times ./ 3600 ./ 24),
    0.0,
    maximum(SWE),
)

lines!(ax2, model_times ./ 3600 ./ 24, SWE, label = "Model", color = "blue")
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
if !comparison_data.P.absent
    lines!(
        ax3,
        data_times ./ 3600 ./ 24,
        (comparison_data.P.values .* (1e3 * 24 * 3600))[data_id_post_spinup],
        label = "Total precip (data)",
    )
end


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

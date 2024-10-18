import SciMLBase
import ClimaTimeSteppers as CTS
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
import ClimaParams as CP
using Plots
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
const FT = Float64
earth_param_set = LP.LandParameters(FT)
climaland_dir = pkgdir(ClimaLand)

include(joinpath(climaland_dir, "experiments/integrated/fluxnet/data_tools.jl"))
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

include(
    joinpath(
        climaland_dir,
        "experiments/integrated/fluxnet/met_drivers_FLUXNET.jl",
    ),
)

# Now we set up the model. For the soil model, we pick
# a model type and model args:
soil_domain = land_domain
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
    PAR_albedo = soil_α_PAR,
    NIR_albedo = soil_α_NIR,
);

soil_args = (domain = soil_domain, parameters = soil_ps)
soil_model_type = Soil.EnergyHydrology{FT}

# Soil microbes model
soilco2_type = Soil.Biogeochemistry.SoilCO2Model{FT}

soilco2_ps = SoilCO2ModelParameters(FT)

Csom = ClimaLand.PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))

# Set the soil CO2 BC to being atmospheric CO2
soilco2_top_bc = Soil.Biogeochemistry.AtmosCO2StateBC()
soilco2_bot_bc = Soil.Biogeochemistry.SoilCO2FluxBC((p, t) -> 0.0) # no flux
soilco2_sources = (MicrobeProduction{FT}(),)

soilco2_boundary_conditions = (; top = soilco2_top_bc, bottom = soilco2_bot_bc)

soilco2_args = (;
    boundary_conditions = soilco2_boundary_conditions,
    sources = soilco2_sources,
    domain = soil_domain,
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
ai_parameterization = PrescribedSiteAreaIndex{FT}(LAIfunction, SAI, RAI)

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

# Integrated plant hydraulics and soil model
land_input = (
    atmos = atmos,
    radiation = radiation,
    soil_organic_carbon = Csom,
    runoff = ClimaLand.Soil.Runoff.SurfaceRunoff(),
)
land = SoilCanopyModel{FT}(;
    soilco2_type = soilco2_type,
    soilco2_args = soilco2_args,
    land_args = land_input,
    soil_model_type = soil_model_type,
    soil_args = soil_args,
    canopy_component_types = canopy_component_types,
    canopy_component_args = canopy_component_args,
    canopy_model_args = canopy_model_args,
)
Y, p, cds = initialize(land)
exp_tendency! = make_exp_tendency(land)
imp_tendency! = make_imp_tendency(land)
jacobian! = make_jacobian(land);
jac_kwargs =
    (; jac_prototype = ClimaLand.ImplicitEquationJacobian(Y), Wfact = jacobian!);

#Initial conditions
Y.soil.ϑ_l =
    drivers.SWC.status != absent ?
    drivers.SWC.values[1 + Int(round(t0 / DATA_DT))] : soil_ν / 2 # Get soil water content at t0
# Both data and simulation are reference to 2005-01-01-00 (LOCAL)
# or 2005-01-01-06 (UTC)
Y.soil.θ_i = FT(0.0)
T_0 =
    drivers.TS.status != absent ?
    drivers.TS.values[1 + Int(round(t0 / DATA_DT))] :
    drivers.TA.values[1 + Int(round(t0 / DATA_DT))] + 40# Get soil temperature at t0
ρc_s =
    volumetric_heat_capacity.(
        Y.soil.ϑ_l,
        Y.soil.θ_i,
        land.soil.parameters.ρc_ds,
        earth_param_set,
    )
Y.soil.ρe_int =
    volumetric_internal_energy.(Y.soil.θ_i, ρc_s, T_0, earth_param_set)
Y.soilco2.C .= FT(0.000412) # set to atmospheric co2, mol co2 per mol air
ψ_stem_0 = FT(-1e5 / 9800) # pressure in the leaf divided by rho_liquid*gravitational acceleration [m]
ψ_leaf_0 = FT(-2e5 / 9800)
ψ_comps = n_stem > 0 ? [ψ_stem_0, ψ_leaf_0] : ψ_leaf_0

S_l_ini =
    inverse_water_retention_curve.(retention_model, ψ_comps, plant_ν, plant_S_s)

for i in 1:(n_stem + n_leaf)
    Y.canopy.hydraulics.ϑ_l.:($i) .=
        augmented_liquid_fraction.(plant_ν, S_l_ini[i])
end

Y.canopy.energy.T = drivers.TA.values[1 + Int(round(t0 / DATA_DT))] # Get atmos temperature at t0

set_initial_cache! = make_set_initial_cache(land)
set_initial_cache!(p, Y, t0);

# Simulation
sv = (;
    t = Array{Float64}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
)
saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat)
## How often we want to update the drivers. Note that this uses the defined `t0` and `tf`
## defined in the simulatons file
updateat = Array(t0:DATA_DT:tf)
model_drivers = ClimaLand.get_drivers(land)
updatefunc = ClimaLand.make_update_drivers(model_drivers)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
cb = SciMLBase.CallbackSet(driver_cb, saving_cb)

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

using ClimaDiagnostics
using ClimaUtilities

outdir = joinpath(pkgdir(ClimaLand), "experiments/integrated/fluxnet/out")
output_dir = ClimaUtilities.OutputPathGenerator.generate_output_path(outdir)

d_writer = ClimaDiagnostics.Writers.DictWriter()

ref_time = DateTime(2005) # random. not sure what it should be

diags = ClimaLand.default_diagnostics(
    land,
    ref_time;
    output_writer = d_writer,
    output_vars = :long,
    average_period = :hourly,
)

diagnostic_handler =
    ClimaDiagnostics.DiagnosticsHandler(diags, Y, p, t0, dt = dt);

diag_cb = ClimaDiagnostics.DiagnosticsCallback(diagnostic_handler);

drivers_diag = ClimaLand.get_drivers(land)
updatefunc = ClimaLand.make_update_drivers(drivers_diag)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)

sol = SciMLBase.solve(
    prob,
    ode_algo;
    dt = dt,
    callback = SciMLBase.CallbackSet(driver_cb, diag_cb),
);

short_names_1D = [
    "sif", # SIF
    "ra", # AR
    "gs", # g_stomata
    "gpp", # GPP
    "ct", # canopy_T
    "swu", # SW_u
    "lwu", # LW_u
    "er", # ER
    "et", # ET
    "msf", # β
    "shf", # SHF
    "lhf", # LHF
    "ghf", # G
    "rn", # Rn
]

short_names_2D = [
    "swc", # swc_sfc or swc_5 or swc_10
    "tsoil", # soil_T_sfc or soil_T_5 or soil_T_10
    "si", # si_sfc or si_5 or si_10
]

hourly_diag_name = short_names_1D .* "_1h_average"
hourly_diag_name_2D = short_names_2D .* "_1h_average"

# diagnostic_as_vectors()[2] is a vector of a variable,
# whereas diagnostic_as_vectors()[1] is a vector or time associated with that variable.
SIF, AR, g_stomata, GPP, canopy_T, SW_u, LW_u, ER, ET, β, SHF, LHF, G, Rn = [
    ClimaLand.Diagnostics.diagnostic_as_vectors(d_writer, diag_name)[2][1:721]
    for diag_name in hourly_diag_name
] # [1:721] only because summer subset of model output

swc_5, soil_T_5, si_5 = [
    ClimaLand.Diagnostics.diagnostic_as_vectors(
        d_writer,
        diag_name;
        layer = 5,
    )[2][1:721] for diag_name in hourly_diag_name_2D
]

times_diag =
    ClimaLand.Diagnostics.diagnostic_as_vectors(d_writer, "swc_1h_average")[1][1:721]

# convert units for GPP and ET
GPP = GPP .* 1e6 # mol to μmol
AR = AR .* 1e6
ET = ET .* 24 .* 3600

# Plotting
daily = sol.t ./ 3600 ./ 24
savedir =
    joinpath(climaland_dir, "experiments/integrated/fluxnet/$site_ID/out/")

if !isdir(savedir)
    mkdir(savedir)
end

# Number of days to plot
num_days = N_days - N_spinup_days

# Time series of model and data outputs
data_times = [0:DATA_DT:(num_days * S_PER_DAY);]
model_times = [0:dt_save:(num_days * S_PER_DAY);]

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
if drivers.GPP.status == absent
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
    GPP_data =
        drivers.GPP.values[Int64(t_spinup ÷ DATA_DT):Int64(tf ÷ DATA_DT)] .* 1e6
    plot_avg_comp(
        "GPP",
        GPP,
        dt_save,
        GPP_data,
        FT(DATA_DT),
        num_days,
        drivers.GPP.units,
        savedir,
    )
end

# Upwelling shortwave radiation is referred to as outgoing in the data
if drivers.SW_OUT.status == absent
    plot_daily_avg("SW up", SW_u, dt_save, num_days, "w/m^2", savedir, "model")
else
    SW_u_data = FT.(drivers.SW_OUT.values)[Int64(t_spinup ÷ DATA_DT):Int64(
        tf ÷ DATA_DT,
    )]
    plot_avg_comp(
        "SW up",
        SW_u,
        dt_save,
        SW_u_data,
        FT(DATA_DT),
        num_days,
        drivers.SW_OUT.units,
        savedir,
    )
end

# Upwelling longwave radiation is referred to outgoing in the data
if drivers.LW_OUT.status == absent
    plot_daily_avg("LW up", LW_u, dt_save, num_days, "w/m^2", savedir, "model")
else
    LW_u_data = FT.(drivers.LW_OUT.values)[Int64(t_spinup ÷ DATA_DT):Int64(
        tf ÷ DATA_DT,
    )]
    plot_avg_comp(
        "LW up",
        LW_u,
        dt_save,
        LW_u_data,
        FT(DATA_DT),
        num_days,
        drivers.LW_OUT.units,
        savedir,
    )
end

# ET
if drivers.LE.status == absent
    plot_daily_avg("ET", ET, dt_save, num_days, "mm/day", savedir, "Model")
else
    measured_T =
        drivers.LE.values ./ (LP.LH_v0(earth_param_set) * 1000) .*
        (1e3 * 24 * 3600)
    ET_data = measured_T[Int64(t_spinup ÷ DATA_DT):Int64(tf ÷ DATA_DT)]
    plot_avg_comp(
        "ET",
        ET,
        dt_save,
        ET_data,
        FT(DATA_DT),
        num_days,
        "mm/day",
        savedir,
    )
end

# Sensible Heat Flux
if drivers.H.status == absent
    plot_daily_avg("SHF", SHF, dt_save, num_days, "w/m^2", savedir, "Model")
else
    SHF_data = drivers.H.values[Int64(t_spinup ÷ DATA_DT):Int64(tf ÷ DATA_DT)]
    plot_avg_comp(
        "SHF",
        SHF,
        dt_save,
        SHF_data,
        FT(DATA_DT),
        N_days - N_spinup_days,
        drivers.H.units,
        savedir,
    )
end

# Latent Heat Flux
if drivers.LE.status == absent
    plot_daily_avg("LHF", LHF, dt_save, num_days, "w/m^2", savedir)
else
    LHF_data = drivers.LE.values[Int64(t_spinup ÷ DATA_DT):Int64(tf ÷ DATA_DT)]
    plot_avg_comp(
        "LHF",
        LHF,
        dt_save,
        LHF_data,
        FT(DATA_DT),
        N_days - N_spinup_days,
        drivers.LE.units,
        savedir,
    )
end

# Ground Heat Flux
# In the land model, positive fluxes are always in the upwards direction
# In the data, a positive G indicates a flux into the soil, so we adjust the sign
# on the data G
if drivers.G.status != absent
    G_data_avg = compute_diurnal_avg(
        FT.(drivers.G.values)[Int64(t_spinup ÷ DATA_DT):Int64(tf ÷ DATA_DT)],
        data_times,
        num_days,
    )
    plt1 = Plots.plot(0.5:0.5:24, -1 .* G_data_avg, label = "Data: -G")
    Plots.plot!(
        plt1,
        ylabel = "Flux (W/m^2)",
        title = "Energy balance at the site",
    )
    if drivers.LE.status != absent &&
       drivers.H.status != absent &&
       drivers.SW_OUT.status != absent &&
       drivers.LW_OUT.status != absent
        Rn_data =
            (drivers.SW_IN.values .- drivers.SW_OUT.values) .+
            (drivers.LW_IN.values .- drivers.LW_OUT.values)
        G_alternate_data_avg = compute_diurnal_avg(
            FT.(drivers.H.values .+ drivers.LE.values .- Rn_data)[Int64(
                t_spinup ÷ DATA_DT,
            ):Int64(tf ÷ DATA_DT)],
            data_times,
            num_days,
        )
        HplusL_avg = compute_diurnal_avg(
            FT.(drivers.H.values .+ drivers.LE.values)[Int64(
                t_spinup ÷ DATA_DT,
            ):Int64(tf ÷ DATA_DT)],
            data_times,
            num_days,
        )
        RminusG_avg = compute_diurnal_avg(
            FT.(Rn_data .- drivers.G.values)[Int64(t_spinup ÷ DATA_DT):Int64(
                tf ÷ DATA_DT,
            )],
            data_times,
            num_days,
        )
        Plots.plot!(
            plt1,
            0.5:0.5:24,
            G_alternate_data_avg,
            label = "Data: (H+L-Rn)_site",
            xlabel = "Hour of day",
        )
        plt2 = Plots.scatter(
            HplusL_avg,
            RminusG_avg,
            label = "Diurnally averaged data",
            xlabel = "H+L",
            ylabel = "R-G",
        )
        Plots.plot(plt1, plt2, layout = (2, 1))
    end
end

Plots.savefig(joinpath(savedir, "energy_balance_data.png"))


# Water stress factor
plot_daily_avg("moisture_stress", β, dt_save, num_days, "", savedir, "Model")

# Stomatal conductance
plot_daily_avg(
    "stomatal conductance",
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

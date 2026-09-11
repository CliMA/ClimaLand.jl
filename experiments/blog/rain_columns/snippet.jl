# Where does the rain go? A 2 m soil column, bare or under a forest,
# hit by a 50 mm storm on the first night of an idealized July, then 29 dry days.
using Dates, ClimaLand, ClimaLand.Soil, ClimaLand.Canopy, CairoMakie
using ClimaLand.Domains: Column, obtain_surface_domain
using ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
import ClimaLand.Parameters as LP, ClimaLand.Simulations as Sim
import ClimaParams, ClimaDiagnostics, Thermodynamics as TD, Insolation

FT = Float64
toml_dict = LP.create_toml_dict(FT)
params = LP.LandParameters(toml_dict)

# --- the place and the weather
long, lat = -92.2, 38.7                 # Missouri: sets the sun's path and default parameter maps
domain = Column(; zlim = (-2.0, 0.0), nelements = 30, dz_tuple = (0.2, 0.025), longlat = (long, lat))
start_date = DateTime(2010, 7, 1, 6)    # 00:00 local time
seconds(t) = float(t)                   # simulation time → seconds since start_date
raining(t) = seconds(t) < 12 * 3600     # 50 mm between midnight and noon on day 1
precip(t) = raining(t) ? -50e-3 / (12 * 3600) : 0.0
T_air(t) = 298.15 + 6 * sin(2π * (seconds(t) / 3600 - 9) / 24)   # 19–31 °C, warmest at 15:00
function q_air(t)                                            # 50 % relative humidity (95 % in the rain)
    e = (raining(t) ? 0.95 : 0.5) * TD.saturation_vapor_pressure(LP.thermodynamic_parameters(params), T_air(t), TD.Liquid())
    return 0.622e / (101325 - 0.378e)
end
cosθ(t) = max(0, Insolation.insolation(start_date + Second(round(Int, seconds(t))), lat, long, params.insol_params).μ)
SW_d(t) = (raining(t) ? 300 : 1000) * cosθ(t)
LW_d(t) = 0.85 * 5.67e-8 * T_air(t)^4
atmos = ClimaLand.PrescribedAtmosphere(TimeVaryingInput(precip), TimeVaryingInput(t -> 0.0),
    TimeVaryingInput(T_air), TimeVaryingInput(t -> 2.0), TimeVaryingInput(q_air),
    TimeVaryingInput(t -> 101325.0), start_date, FT(20), toml_dict)
radiation = ClimaLand.PrescribedRadiativeFluxes(FT, TimeVaryingInput(SW_d), TimeVaryingInput(LW_d), start_date;
    toml_dict, cosθs = (t, s) -> ClimaLand.default_cos_zenith_angle(t, s; insol_params = params.insol_params, longitude = long, latitude = lat))
forcing = (; atmos, radiation)

# --- the soil: a loam (van Genuchten parameters from Carsel & Parrish, 1988)
loam = (; ν = 0.43, θ_r = 0.078, K_sat = 2.89e-6, hydrology_cm = vanGenuchten{FT}(; α = 3.6, n = 1.56))
soil(; kw...) = Soil.EnergyHydrology{FT}(domain, forcing, toml_dict; retention_parameters = loam,
    composition_parameters = (; ν_ss_om = 0.0, ν_ss_quartz = 0.4, ν_ss_gravel = 0.0), S_s = 1e-3,
    runoff = Soil.Runoff.SurfaceRunoff(), bottom_bc = Soil.EnergyWaterFreeDrainage(), kw...)
bare = soil()

# --- the same soil under a forest: 5 m² of leaves per m² of ground, roots reaching ~1 m
components = (:canopy, :soil, :soilco2)
surface = obtain_surface_domain(domain)
LAI = TimeVaryingInput(t -> 5.0)
canopy = Canopy.CanopyModel{FT}(surface, (; atmos, radiation, ground = ClimaLand.PrognosticGroundConditions{FT}()), LAI, toml_dict;
    prognostic_land_components = components,
    biomass = Canopy.PrescribedBiomassModel{FT}(surface, LAI, toml_dict; rooting_depth = 1.0, height = 2.0),
    hydraulics = Canopy.PlantHydraulicsModel{FT}(surface, toml_dict; retention_model = Canopy.LinearRetentionCurve{FT}(5e-5)),
    photosynthesis = Canopy.FarquharModel{FT}(surface, toml_dict),
    conductance = Canopy.MedlynConductanceModel{FT}(surface, toml_dict))
forest = SoilCanopyModel{FT}(forcing, LAI, toml_dict, domain;
    soil = soil(; prognostic_land_components = components, additional_sources = (ClimaLand.RootExtraction{FT}(),)), canopy)

# --- initial state: every layer at a suction of 2 m (≈ −20 kPa), the plant in equilibrium with it
function set_ic!(Y, p, t0, model)
    s = model isa SoilCanopyModel ? model.soil : model
    (; ν, θ_r, hydrology_cm, ρc_ds) = s.parameters
    Y.soil.ϑ_l .= θ_r + Soil.inverse_matric_potential(hydrology_cm, -2.0) * (ν - θ_r)
    Y.soil.θ_i .= 0
    ρc = Soil.volumetric_heat_capacity.(Y.soil.ϑ_l, Y.soil.θ_i, ρc_ds, params)
    Y.soil.ρe_int .= Soil.volumetric_internal_energy.(Y.soil.θ_i, ρc, T_air(0.0), params)
    if model isa SoilCanopyModel
        for (dst, src) in ((p.drivers.T, atmos.T), (p.drivers.P, atmos.P), (p.drivers.q, atmos.q), (p.drivers.c_co2, atmos.c_co2))
            ClimaUtilities.TimeVaryingInputs.evaluate!(dst, src, t0)
        end
        Sim.set_soilco2_initial_conditions!(Y, p, model)
        h = model.canopy.hydraulics.parameters
        ψ_plant = -2.0 - (model.canopy.biomass.height / 2 + model.canopy.biomass.rooting_depth)
        Y.canopy.hydraulics.ϑ_l .= Canopy.augmented_liquid_fraction(h.ν, Canopy.inverse_water_retention_curve(h.retention_model, ψ_plant, h.ν, h.S_s))
        Y.canopy.energy.T .= T_air(0.0)
    end
end
import ClimaUtilities

# --- run both for 30 days, saving hourly means
function run(model, vars)
    writer = ClimaDiagnostics.Writers.DictWriter()
    diagnostics = ClimaLand.default_diagnostics(model, start_date; output_writer = writer, output_vars = vars, reduction_period = :hourly)
    sim = Sim.LandSimulation(start_date, start_date + Day(30), 60.0, model; set_ic!, updateat = Second(60), user_callbacks = (), diagnostics)
    Sim.solve!(sim)
    return writer
end
out_bare = run(bare, ["swc", "et"])
out_forest = run(forest, ["swc", "et"])

# --- plot: the bare column through time, both columns' profiles on day 30, and the water returned to the air
z = vec(parent(domain.fields.z))
swc(w) = (d = w["swc_1h_average"]; ts = sort(collect(keys(d))); (float.(ts) ./ 86400, reduce(hcat, vec(parent(d[t])) for t in ts)))
et(w, scale) = (d = w["et_1h_average"]; ts = sort(collect(keys(d))); (float.(ts) ./ 86400, [only(parent(d[t])) * scale for t in ts]))
days, θ_bare = swc(out_bare); _, θ_forest = swc(out_forest)
θ0 = loam.θ_r + Soil.inverse_matric_potential(loam.hydrology_cm, -2.0) * (loam.ν - loam.θ_r)
fig = Figure(size = (1000, 700))
ax = Axis(fig[1, 1]; title = "bare loam", xlabel = "days since the storm", ylabel = "depth (m)")
hm = heatmap!(ax, days, z, permutedims(θ_bare); colormap = :YlGnBu, colorrange = (0, 0.45))
Colorbar(fig[1, 2], hm; label = "soil water content (m³/m³)")
ax = Axis(fig[1, 3]; title = "profiles on day 30", xlabel = "soil water content (m³/m³)")
vlines!(ax, [θ0]; color = :gray, linestyle = :dash, label = "before the storm")
lines!(ax, θ_bare[:, end], z; linewidth = 3, label = "bare loam")
lines!(ax, θ_forest[:, end], z; linewidth = 3, label = "loam + forest")
axislegend(ax; position = :lb)
ax = Axis(fig[2, 1:3]; xlabel = "days since the storm", ylabel = "water to the air (mm/day)")
lines!(ax, et(out_bare, 8.64e7)...; label = "bare loam", linewidth = 2)       # soil model: m/s → mm/day
lines!(ax, et(out_forest, 86400.0)...; label = "loam + forest", linewidth = 2) # coupled model: kg/m²/s → mm/day
axislegend(ax)
save("where_does_the_rain_go.png", fig)

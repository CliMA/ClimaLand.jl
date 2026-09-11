# Where does the rain go? One storm on five columns: sand, loam, clay (bare),
# and loam under grass and under forest. Idealized July weather, 30-day dry-down.
using Dates, Statistics, Serialization
import ClimaParams as CP
import ClimaDiagnostics
import Insolation
using ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
using ClimaLand
using ClimaLand.Domains: Column, obtain_surface_domain
using ClimaLand.Soil, ClimaLand.Canopy
import ClimaLand.Simulations: LandSimulation, solve!
import ClimaLand.Parameters as LP
import Thermodynamics as TD

const FT = Float64
toml_dict = LP.create_toml_dict(FT)
earth_param_set = LP.LandParameters(toml_dict)
thermo_params = LP.thermodynamic_parameters(earth_param_set)

NDAYS = parse(Int, get(ENV, "NDAYS", "30"))
CASES = split(get(ENV, "CASES", "sand_bare,loam_bare,clay_bare,loam_grass,loam_forest"), ",")
OUTDIR = get(ENV, "OUTDIR", "out")
DT = parse(Float64, get(ENV, "DT", "60"))
PLANT_A = parse(Float64, get(ENV, "PLANT_A", "5e-5"))
mkpath(OUTDIR)

# ---------------------------------------------------------------- domain
long, lat = FT(-92.2), FT(38.7)              # Missouri: sets sun angle & default maps
domain = Column(; zlim = (FT(-2), FT(0)), nelements = 30, dz_tuple = (FT(0.2), FT(0.025)), longlat = (long, lat))
surface_domain = obtain_surface_domain(domain)

# ---------------------------------------------------------------- idealized forcing
start_date = DateTime(2010, 7, 1, 6)          # 00:00 local (UTC-6)
stop_date = start_date + Day(NDAYS)
hour_of_day(t) = mod(float(t) / 3600, 24)
STORM_HOURS = 12                              # 50 mm of rain, midnight to noon on day 1
raining(t) = float(t) < STORM_HOURS * 3600
precip(t) = raining(t) ? -50e-3 / (STORM_HOURS * 3600) : 0.0  # m/s, negative = downward
T_air(t) = 298.15 + 6 * sin(2π * (hour_of_day(t) - 9) / 24)  # 22–34 °C, peak 15:00
RH(t) = raining(t) ? 0.95 : 0.5
function q_air(t)
    e_sat = TD.saturation_vapor_pressure(thermo_params, T_air(t), TD.Liquid())
    e = RH(t) * e_sat
    return 0.622 * e / (101325 - 0.378 * e)
end
cosθ(t) = max(0, Insolation.insolation(start_date + Second(round(Int, float(t))), lat, long, earth_param_set.insol_params).μ)
SW_d(t) = (raining(t) ? 300 : 1000) * cosθ(t)   # overcast while it rains
LW_d(t) = 0.85 * 5.67e-8 * T_air(t)^4

atmos = ClimaLand.PrescribedAtmosphere(
    TimeVaryingInput(precip), TimeVaryingInput(t -> 0.0),  # rain, snow
    TimeVaryingInput(T_air), TimeVaryingInput(t -> 2.0),   # temperature, wind
    TimeVaryingInput(q_air), TimeVaryingInput(t -> 101325.0),
    start_date, FT(20), toml_dict)                        # reference height 20 m
radiation = ClimaLand.PrescribedRadiativeFluxes(FT, TimeVaryingInput(SW_d), TimeVaryingInput(LW_d), start_date; toml_dict,
    cosθs = (t, s) -> ClimaLand.default_cos_zenith_angle(t, s; insol_params = earth_param_set.insol_params, longitude = long, latitude = lat))
forcing = (; atmos, radiation)

# ---------------------------------------------------------------- soils (Carsel & Parrish 1988)
soils = Dict(
    "sand" => (; ν = 0.43, θ_r = 0.045, α = 14.5, n = 2.68, K_sat = 8.25e-5, quartz = 0.9),
    "loam" => (; ν = 0.43, θ_r = 0.078, α = 3.6, n = 1.56, K_sat = 2.89e-6, quartz = 0.4),
    "clay" => (; ν = 0.38, θ_r = 0.068, α = 0.8, n = 1.09, K_sat = 5.56e-7, quartz = 0.2),
)
plants = Dict(
    "grass" => (; LAI = 2.0, height = 0.5, rooting_depth = 0.3),
    "forest" => (; LAI = 5.0, height = 2.0, rooting_depth = 1.0),
)

function make_soil(s; kwargs...)
    retention_parameters = (; ν = FT(s.ν), θ_r = FT(s.θ_r), K_sat = FT(s.K_sat),
        hydrology_cm = vanGenuchten{FT}(; α = FT(s.α), n = FT(s.n)))
    composition_parameters = (; ν_ss_om = FT(0), ν_ss_quartz = FT(s.quartz), ν_ss_gravel = FT(0))
    Soil.EnergyHydrology{FT}(domain, forcing, toml_dict; retention_parameters, composition_parameters,
        S_s = FT(1e-3), runoff = Soil.Runoff.SurfaceRunoff(), bottom_bc = Soil.EnergyWaterFreeDrainage(), kwargs...)
end

function make_model(soilname, cover)
    s = soils[soilname]
    cover == "bare" && return make_soil(s)
    v = plants[cover]
    components = (:canopy, :soil, :soilco2)
    soil = make_soil(s; prognostic_land_components = components, additional_sources = (ClimaLand.RootExtraction{FT}(),))
    LAI = TimeVaryingInput(t -> FT(v.LAI))
    biomass = Canopy.PrescribedBiomassModel{FT}(surface_domain, LAI, toml_dict; rooting_depth = FT(v.rooting_depth), height = FT(v.height))
    # plant water capacitance: the default slope (0.00196 /m) stores ~400 mm per MPa for this forest; real trees hold ~5-20 mm per MPa
    hydraulics = Canopy.PlantHydraulicsModel{FT}(surface_domain, toml_dict; retention_model = Canopy.LinearRetentionCurve{FT}(FT(PLANT_A)))
    canopy = Canopy.CanopyModel{FT}(surface_domain, (; atmos, radiation, ground = ClimaLand.PrognosticGroundConditions{FT}()), LAI, toml_dict;
        prognostic_land_components = components, biomass, hydraulics,
        # Farquhar + Medlyn: stomata respond to leaf water potential within the hour
        # (the default P-model applies water stress through a ~15-day acclimation)
        photosynthesis = Canopy.FarquharModel{FT}(surface_domain, toml_dict),
        conductance = Canopy.MedlynConductanceModel{FT}(surface_domain, toml_dict))
    return SoilCanopyModel{FT}(forcing, LAI, toml_dict, domain; soil, canopy)
end

# ---------------------------------------------------------------- initial conditions: same suction everywhere
ψ0 = FT(-2)  # m of head, ≈ -20 kPa
function set_ic!(Y, p, t0, model)
    soil = model isa Soil.EnergyHydrology ? model : model.soil
    (; ν, θ_r, hydrology_cm, S_s, ρc_ds) = soil.parameters
    S = Soil.inverse_matric_potential(hydrology_cm, ψ0)
    Y.soil.ϑ_l .= θ_r + S * (ν - θ_r)
    Y.soil.θ_i .= 0
    T0 = FT(T_air(0.0))
    ρc_s = Soil.volumetric_heat_capacity.(Y.soil.ϑ_l, Y.soil.θ_i, ρc_ds, earth_param_set)
    Y.soil.ρe_int .= Soil.volumetric_internal_energy.(Y.soil.θ_i, ρc_s, T0, earth_param_set)
    if model isa SoilCanopyModel
        for (dst, src) in ((p.drivers.T, atmos.T), (p.drivers.P, atmos.P), (p.drivers.q, atmos.q), (p.drivers.c_co2, atmos.c_co2))
            ClimaUtilities.TimeVaryingInputs.evaluate!(dst, src, t0)
        end
        ClimaLand.Simulations.set_soilco2_initial_conditions!(Y, p, model)
        hyd = model.canopy.hydraulics.parameters
        # plant in hydraulic equilibrium with the soil: soil suction plus the lift from roots to leaves
        ψ_plant0 = ψ0 - (model.canopy.biomass.height / 2 + model.canopy.biomass.rooting_depth)
        S_l = Canopy.inverse_water_retention_curve(hyd.retention_model, ψ_plant0, hyd.ν, hyd.S_s)
        Y.canopy.hydraulics.ϑ_l .= Canopy.augmented_liquid_fraction(hyd.ν, S_l)
        Y.canopy.energy.T .= T0
        model.canopy.photosynthesis isa Canopy.PModel &&
            ClimaLand.Simulations.set_canopy_component_initial_conditions!(Y, p, model.canopy.photosynthesis, model.canopy)
    end
end
import ClimaUtilities

# ---------------------------------------------------------------- run
function run_case(name)
    soilname, cover = split(name, "_")
    model = make_model(String(soilname), String(cover))
    possible = ClimaLand.Diagnostics.get_possible_diagnostics(model)
    wanted = ["swc", "swp", "tsoil", "infil", "sr", "sdr", "et", "lhf", "shf", "trans", "soillhf", "lwp", "msf", "gpp", "rn", "far"]
    output_vars = filter(v -> v in possible, wanted)
    writer = ClimaDiagnostics.Writers.DictWriter()
    diagnostics = ClimaLand.default_diagnostics(model, start_date; output_writer = writer, output_vars, reduction_period = :hourly)
    sim = LandSimulation(start_date, stop_date, DT, model; set_ic!, updateat = Second(round(Int, DT)), user_callbacks = (), diagnostics)
    t = @elapsed solve!(sim)
    @info "$name done in $(round(t; digits = 1)) s"
    out = Dict{String, Any}("z" => vec(parent(domain.fields.z)))
    for (k, v) in writer.dict
        times = sort(collect(keys(v)))
        out[k] = (; times = [float(τ) for τ in times], data = [vec(parent(v[τ])) for τ in times])
    end
    return out
end

results = Dict(name => run_case(name) for name in CASES)
serialize(joinpath(OUTDIR, "results.jls"), results)
println("saved ", joinpath(OUTDIR, "results.jls"))

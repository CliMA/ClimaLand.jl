"""
test_minimal_forward.jl — minimaler Forward-Lauf zum Eingrenzen des
`DomainError: log ... negative real argument` (aus `q_vap_saturation`), der seit
dem Rebase auf das neue main auftritt.

STAND DER DIAGNOSE (Stufe 0 = reiner Default-Aufbau):
    HARV, LENO, WOOD   PASS
    MOAB, JORN         FAIL   (trocken, LAI 0.24-0.32)
    LENO Stufe 2       FAIL   (Moisture-Stress zugeschaltet)
Also ZWEI Effekte: trockene Standorte kippen schon im Default, feuchte erst mit
PiecewiseMoistureStressModel.

HAUPTVERDACHT (PR #1774, canopy_boundary_fluxes.jl): `q_canopy` wird jetzt VOR
dem Monin-Obukhov-Solve extern gebildet statt im Solver
    -  q_canopy = inputs.q_vap_sfc_guess
    +  q_canopy = component_specific_humidity(model, Y, p)   # = q_sat(T_canopy)
`component_specific_humidity` ruft genau die Funktion aus dem Stacktrace auf.
Dieser Pfad ist auf JEDER Stufe aktiv, koennte also beide Effekte erklaeren.

Dieses Skript schreibt HALBSTUENDLICHE Diagnostics in eine CSV (auch wenn der
Lauf crasht) und prueft zusaetzlich per Callback jeden Schritt auf nicht-endliche
oder unphysikalische Temperaturen, damit die Absturzstelle sichtbar wird. Die
Tagesmittel der Pipeline verstecken den Weglauf, weil er in wenigen Schritten
passiert.

AUFRUF:
    julia --project=.buildkite experiments/calibrate_neon_pipeline_AllDepth/test_minimal_forward.jl

    SITE=NEON-jorn LEVEL=0 YEARS=1 julia --project=.buildkite \\
        experiments/calibrate_neon_pipeline_AllDepth/test_minimal_forward.jl

ENV-VARIABLEN: SITE (NEON-harv), LEVEL (0), YEARS (1), DT (900.0), START (2018-01-01)

STUFEN (LEVEL):
    0  alles Default  (== experiments/calibration/models/snowy_land.jl)
    1  + PModel + PModelConductance
    2  + PiecewiseMoistureStressModel, mit Gupta-retention_parameters, die
       Stress-Funktion UND Bodenmodell teilen (LandModel verlangt das:
       check_land_equality, land.jl:95-102, exaktes == auf θ_high/ν)
    3  wie 2, aber Rosetta statt Gupta als Bodenkarte

CSV: /kiwi-data/Data/groupMembers/evametz/ClimaLand_Output/Test/PR1774/
     diag_<site>_L<level>_dt<dt>.csv
"""

import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaLand.Simulations: LandSimulation, solve!
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Canopy
using ClimaLand.Snow
using ClimaCore
using ClimaDiagnostics
import ClimaUtilities.TimeManager: date, ITime
using Dates
using DataFrames
using CSV
using Printf

include(joinpath(@__DIR__, "src/site_metadata.jl"))

const FT = Float64
const OUTDIR = "/kiwi-data/Data/groupMembers/evametz/ClimaLand_Output/Test/PR1774"

site  = get(ENV, "SITE", "NEON-harv")
level = parse(Int, get(ENV, "LEVEL", "0"))
years = parse(Int, get(ENV, "YEARS", "1"))
DT    = FT(parse(Float64, get(ENV, "DT", "900.0")))

start_date = DateTime(get(ENV, "START", "2018-01-01"))
stop_date  = start_date + Year(years)

mkpath(OUTDIR)
println("="^72)
println("site=$site  level=$level  dt=$(DT)s  $(Date(start_date)) -> $(Date(stop_date))")
println("out=$OUTDIR")
println("="^72)

md = _get_neon_site_metadata(site)
lat, long = FT(md.lat), FT(md.long)

# Domain identisch zu ForwardRun.jl
domain = Column(; zlim = (FT(-6.2), FT(0)), nelements = 24,
    dz_tuple = (FT(2), FT(0.038)), longlat = (long, lat))
canopy_domain = ClimaLand.Domains.obtain_surface_domain(domain)
surface_space = domain.space.surface

toml_dict = LP.create_toml_dict(FT)

(; atmos, radiation) = FluxnetSimulations.prescribed_forcing_fluxnet(
    site, lat, long, 0, FT(md.atmos_h), start_date, toml_dict, FT)
forcing = (; atmos, radiation)

LAI = ClimaLand.Canopy.prescribed_climatological_lai_modis(surface_space)
prognostic_land_components = (:canopy, :snow, :soil, :soilco2)

# ── Komponenten je nach Stufe ────────────────────────────────────────────────
canopy_kwargs = Dict{Symbol, Any}()
if level >= 1
    println("  + PModel / PModelConductance")
    canopy_kwargs[:photosynthesis] = PModel{FT}(domain, toml_dict)
    canopy_kwargs[:conductance] = PModelConductance{FT}(toml_dict)
end

# retention_parameters wird von Stufe 2 UND 3 gebraucht und MUSS dasselbe Objekt
# sein: LandModels check_land_equality vergleicht θ_high == ν mit exaktem ==.
retention_parameters = nothing
if level >= 2
    # Stufe 2: Gupta, Stufe 3: Rosetta — jeweils EIN Objekt, das Stress-Funktion
    # und Bodenmodell teilen (check_land_equality verlangt exaktes ==).
    src = level >= 3 ? "Rosetta" : "Gupta"
    println("  + PiecewiseMoistureStressModel (geteilte $src-retention_parameters)")
    retention_parameters = level >= 3 ?
        Soil.rosetta_soil_vangenuchten_parameters(domain.space.subsurface, FT) :
        Soil.soil_vangenuchten_parameters(domain.space.subsurface, FT)
    canopy_kwargs[:soil_moisture_stress] =
        ClimaLand.Canopy.PiecewiseMoistureStressModel{FT}(
            domain, toml_dict; soil_params = retention_parameters)
end

canopy = ClimaLand.Canopy.CanopyModel{FT}(
    canopy_domain,
    (; atmos, radiation, ground = ClimaLand.PrognosticGroundConditions{FT}()),
    LAI, toml_dict; prognostic_land_components, canopy_kwargs...)

snow = Snow.SnowModel(FT, canopy_domain, forcing, toml_dict, DT;
    prognostic_land_components,
    α_snow = Snow.ZenithAngleAlbedoModel(toml_dict))

soil_kwargs = Dict{Symbol, Any}(
    :prognostic_land_components => prognostic_land_components,
    :additional_sources => (ClimaLand.RootExtraction{FT}(),))
# Ab Stufe 2 MUSS das Bodenmodell dasselbe retention_parameters-Objekt bekommen
# wie die Stress-Funktion: LandModel prueft in check_land_equality (land.jl:95-102)
# θ_high == soil.parameters.ν und θ_low == soil.parameters.θ_r mit exaktem ==.
# Sonst -> AssertionError: tmp1 == tmp2.
level >= 2 && (soil_kwargs[:retention_parameters] = retention_parameters)
soil = Soil.EnergyHydrology{FT}(domain, forcing, toml_dict; soil_kwargs...)

land = LandModel{FT}(forcing, LAI, toml_dict, domain, DT;
    prognostic_land_components, snow, canopy, soil)

# ── Diagnostics: halbstuendlich ─────────────────────────────────────────────
# NICHT dabei: `rsc` (existiert nur auf dem _canopyissue-Branch) und `an`
# (in ClimaLand registriert, aber compute_photosynthesis_net_leaf! ist
# upstream NICHT definiert -> UndefVarError). `gpp` deckt die Photosynthese ab.
output_vars = [
    # A: die drei q_vap_saturation-Nutzer + Referenz
    "ct", "tair", "tsoil", "snowtsfc",
    # B: Pfad zur Canopy-Feuchte
    "msf", "gpp", "lai", "swc",
    # C: Canopy-Energiebilanz
    "clhf", "cshf", "swn", "lwn",
    # D: Wasser + Schnee
    "trans", "et", "swe", "snowc",
]
writer = ClimaDiagnostics.Writers.DictWriter()
diags = ClimaLand.default_diagnostics(land, start_date;
    output_writer = writer, output_vars, reduction_period = :halfhourly)

set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
    site, start_date, 0, land)

# ── Callback: jeden Schritt auf Weglauf pruefen ─────────────────────────────
# Der Absturz passiert INNERHALB eines Zeitschritts, bevor Diagnostics greifen.
# Dieser Callback laeuft nach jedem Schritt und meldet den ERSTEN Schritt, in dem
# eine prognostische Variable nicht-endlich wird oder T unphysikalisch ist.
first_bad = Ref{Union{Nothing, NamedTuple}}(nothing)
function check_state!(integrator)
    first_bad[] === nothing || return nothing
    Y = integrator.u
    for grp in propertynames(Y)
        sub = getproperty(Y, grp)
        for v in propertynames(sub)
            arr = parent(getproperty(sub, v))
            if !all(isfinite, arr)
                first_bad[] = (; t = float(integrator.t), var = "$grp.$v",
                               reason = "nicht endlich",
                               vmin = NaN, vmax = NaN)
                return nothing
            end
            # Temperaturen: negative absolute T ist der Crash-Auslöser
            if v in (:T,) && minimum(arr) <= 0
                first_bad[] = (; t = float(integrator.t), var = "$grp.$v",
                               reason = "T <= 0 K",
                               vmin = minimum(arr), vmax = maximum(arr))
                return nothing
            end
        end
    end
    return nothing
end
# IntervalBasedCallback(period, t0, dt, affect!) — genau wie der eingebaute
# NaNCheckCallback (shared_utilities/utils.jl:684). t0/Δt als ITime, wie es die
# LandSimulation-Konstruktoren selbst tun (Simulations.jl:296-303).
t0_itime = ITime(0, epoch = start_date)
dt_itime = ITime(round(Int, DT), epoch = start_date)
check_cb = ClimaLand.IntervalBasedCallback(
    Second(round(Int, DT)), t0_itime, dt_itime, check_state!)

simulation = LandSimulation(start_date, stop_date, DT, land;
    set_ic! = set_ic!, updateat = Second(DT), diagnostics = diags,
    user_callbacks = (check_cb,))

println("\nlaeuft ...")
@time solve!(simulation)

# ── CSV schreiben (auch nach Crash) ─────────────────────────────────────────
# ACHTUNG Layer-Indizierung: in diagnostic_as_vectors (diagnostics/diagnostic.jl:204)
# ist Index 1 die UNTERSTE Schicht (z = -6.2 m), nlayers(field) die oberste.
# `layer = nothing` liefert die oberste Schicht — genau was wir fuer tsoil/swc
# wollen. Oberflaechenvariablen haben nur eine Schicht, dort ist es ohnehin egal.
# Fehler werden GEMELDET, nicht verschluckt, damit ein fehlender Name auffaellt.
function series(name)
    try
        (t, v) = ClimaLand.Diagnostics.diagnostic_as_vectors(writer, name)
        d = t isa Vector{DateTime} ? t : date.(t)
        return DataFrame(datetime = d, value = Float64.(v))
    catch e
        @warn "konnte $name nicht lesen" exception = e
        return nothing
    end
end

df = nothing
for v in output_vars
    global df                     # sonst behandelt Julia `df` in der Schleife als lokal
    s = series("$(v)_30m_average")
    s === nothing && (@warn "keine Daten fuer $v"; continue)
    rename!(s, :value => Symbol(v))
    df = df === nothing ? s : outerjoin(df, s, on = :datetime)
end

csv_path = joinpath(OUTDIR,
    @sprintf("diag_%s_L%d_dt%d.csv", replace(site, "NEON-" => ""), level, Int(DT)))
if df !== nothing
    sort!(df, :datetime)
    CSV.write(csv_path, df)
    println("\nCSV: $csv_path  ($(nrow(df)) Zeilen, $(ncol(df)) Spalten)")
else
    println("\nWARNUNG: keine Diagnostics erhalten, keine CSV geschrieben")
end

# ── Ergebnis ────────────────────────────────────────────────────────────────
integ = simulation._integrator
t_end = float(integ.t)
t_soll = float(integ.sol.prob.tspan[2])
reached = t_end >= t_soll - float(DT)
finite = all(all(isfinite, parent(getproperty(integ.u, k)))
             for k in propertynames(integ.u))

println()
println("="^72)
if first_bad[] !== nothing
    b = first_bad[]
    println("  ERSTER WEGLAUF:")
    println("    t        = $(b.t) s  (= $(start_date + Second(round(Int, b.t))))")
    println("    Variable = $(b.var)")
    println("    Grund    = $(b.reason)")
    isnan(b.vmin) || println("    min/max  = $(b.vmin) / $(b.vmax)")
else
    println("  kein Weglauf im Callback registriert")
end
println("  Endzeit erreicht : $reached   (t=$t_end von $t_soll)")
println("  Zustand endlich  : $finite")
println("  ERGEBNIS         : ", (reached && finite) ? "PASS" : "FAIL")
println("="^72)

exit((reached && finite) ? 0 : 1)

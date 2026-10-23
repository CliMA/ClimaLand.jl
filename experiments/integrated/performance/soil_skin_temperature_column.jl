# # Soil skin temperature at a jump of the ground conductance

# A single column at (14°N, 8°W), in the Sahel, where the global P-model long
# run of #1907 produced NaNs (#1922). The column uses the default `LandModel`
# with MODIS LAI, as that run does, so its canopy is 0.5 m tall (CLM) with stem
# and root area indices from the MODIS maximum LAI, and the parameter values of
# that run, from before the recalibration of #1926
# (`soil_skin_temperature_column.toml`). An idealized forcing recreates the
# conditions of the failure: a 6 mm/h shower under full sun wets a hot soil
# below canopy air, at a sub-canopy wind of 0.13 m/s.

# There, the heat conductance of the ground jumps between its capped stable
# value and a much larger unstable value where the virtual temperatures of the
# soil surface and the canopy air cross. The Newton updates of the skin
# temperature within the Monin-Obukhov iterations neglect the dependence of the
# conductance on the skin temperature, and can end on an iterate far from the
# skin energy balance; `Soil.solve_soil_surface_temperature_at_a_point` then
# refines the skin temperature by bisection. The column is run with and without
# that refinement; the latter raises `Soil.SKIN_ENERGY_BALANCE_TOLERANCE`, which
# requires Julia 1.12 or later. Without it, at the end of the shower, the soil
# latent heat flux reaches ~25 kW/m², dries the top cell below its residual
# water content, and the run fails. With it, the column runs through.

# Simulation Setup
# Domain: column, 15 m deep, 15 layers from 3 m thick at depth to 5 cm at the top
# Simulation duration: 2 days
# Timestep: 900 s
# Timestepper: ARS111
# Atmospheric and radiation driver updates: Every 3 hours

import ClimaComms
ClimaComms.@import_required_backends
using CairoMakie
using ClimaCore
using Dates
import Insolation
import SurfaceFluxes
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
import ClimaUtilities.OutputPathGenerator: generate_output_path
using ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.Simulations: LandSimulation, step!

const FT = Float64
start_date = DateTime(2009, 10, 15)
stop_date = start_date + Day(2)
Δt = 900.0
longlat = FT.((-8, 14))
shower_rate = 6.0 # mm/h
outdir = generate_output_path(
    "experiments/integrated/performance/soil_skin_temperature_column",
)
toml_dict = LP.create_toml_dict(
    FT;
    override_files = [joinpath(@__DIR__, "soil_skin_temperature_column.toml")],
)

"""
    sahel_column(toml_dict)

Return the default `LandModel` of the column, forced by a hot day with a
shower of `shower_rate` from 12:00 to 15:00 UTC on the first day, under full
sun.
"""
function sahel_column(toml_dict)
    domain = ClimaLand.Domains.Column(;
        zlim = FT.((-15, 0)),
        longlat,
        nelements = 15,
        dz_tuple = FT.((3, 0.05)),
    )
    insol_params = LP.LandParameters(toml_dict).insol_params
    hour(t) = mod(float(t) / 3600, 24)
    raining(t) = float(t) < 86400 && 12 <= hour(t) < 15
    cos_zenith(t) = Insolation.insolation(
        start_date + Second(round(Int, float(t))),
        longlat[2],
        longlat[1],
        insol_params,
    ).μ
    atmos = ClimaLand.PrescribedAtmosphere(
        TimeVaryingInput(t -> raining(t) ? -shower_rate / 3.6e6 : 0.0),
        TimeVaryingInput(t -> 0.0),
        TimeVaryingInput(t -> 304 + 7 * sin(2π * (hour(t) - 9.5) / 24)),
        TimeVaryingInput(t -> 1.0), # the minimum wind speed of the gustiness
        TimeVaryingInput(t -> 0.017),
        TimeVaryingInput(t -> 97000.0),
        start_date,
        FT(10),
        toml_dict;
        reference = SurfaceFluxes.ReferenceAboveApparentSink(),
    )
    radiation = ClimaLand.PrescribedRadiativeFluxes(
        FT,
        TimeVaryingInput(t -> 1000 * max(cos_zenith(t), 0)),
        TimeVaryingInput(t -> 420.0),
        start_date;
        cosθs = (t, start_date) -> cos_zenith(t),
        frac_diff = TimeVaryingInput(t -> 0.3),
        toml_dict,
    )
    LAI = ClimaLand.Canopy.prescribed_lai_modis(
        domain.space.surface,
        start_date,
        stop_date,
    )
    return LandModel{FT}((; atmos, radiation), LAI, toml_dict, domain, Δt)
end

root_cause(err::TaskFailedException) = root_cause(err.task.exception)
root_cause(err::CompositeException) = root_cause(first(err.exceptions))
root_cause(err) = err

"""
    run_column(toml_dict)

Run the column and return, at each step, the time [hours], the soil latent heat
flux [W/m²], and the water content of the top cell above its residual value,
with the error or NaN that stopped the run, if any.
"""
function run_column(toml_dict)
    model = sahel_column(toml_dict)
    simulation = LandSimulation(
        start_date,
        stop_date,
        Δt,
        model;
        outdir,
        diagnostics = nothing,
        user_callbacks = (),
    )
    Y = simulation._integrator.u
    p = simulation._integrator.p
    θ_r_top = Array(parent(model.soil.parameters.θ_r))[end]
    hours = FT[]
    lhf = FT[]
    θ_above_residual = FT[]
    failure = nothing
    for n in 1:round(Int, Dates.value(Second(stop_date - start_date)) / Δt)
        try
            step!(simulation)
        catch err
            failure = root_cause(err)
            break
        end
        if any(isnan, parent(Y.soil.ϑ_l)) || any(isnan, parent(Y.soil.ρe_int))
            failure = ErrorException("NaN in the soil state")
            break
        end
        push!(hours, n * Δt / 3600)
        push!(lhf, only(Array(parent(p.soil.turbulent_fluxes.lhf))))
        push!(θ_above_residual, Array(parent(Y.soil.ϑ_l))[end] - θ_r_top)
    end
    return (; hours, lhf, θ_above_residual, failure)
end

function report(label, run)
    i = argmax(run.lhf)
    max_soil_lhf = run.lhf[i]
    at_hour = run.hours[i]
    min_top_cell_water_above_residual = minimum(run.θ_above_residual)
    last_hour = run.hours[end]
    failure =
        isnothing(run.failure) ? "none" :
        first(split(sprint(showerror, run.failure), '\n'))
    @info label max_soil_lhf at_hour min_top_cell_water_above_residual last_hour failure
end

with_bisection = run_column(toml_dict)
report("With the bisection", with_bisection)

# Without the bisection: the energy balance tolerance is never exceeded
if VERSION >= v"1.12"
    @eval ClimaLand.Soil const SKIN_ENERGY_BALANCE_TOLERANCE = Inf
    without_bisection = Base.invokelatest(run_column, toml_dict)
    report("Without the bisection", without_bisection)
else
    @info "Skipping the run without the bisection, which requires Julia 1.12"
    without_bisection = nothing
end

fig = Figure(size = (800, 600))
ax_lhf = Axis(
    fig[1, 1],
    ylabel = "Soil latent heat flux (W/m²)",
    title = "Sahel column (14°N, 8°W); shower from 12 to 15 UTC on day 1",
)
ax_θ = Axis(
    fig[2, 1],
    xlabel = "Hours since $start_date",
    ylabel = "Top-cell water above residual (m³/m³)",
)
for ax in (ax_lhf, ax_θ)
    vspan!(ax, 12, 15, color = (:steelblue, 0.15))
end
hlines!(ax_θ, 0, color = :black, linestyle = :dash)
for (label, run) in (
    ("with bisection", with_bisection),
    ("without bisection", without_bisection),
)
    isnothing(run) && continue
    lines!(ax_lhf, run.hours, run.lhf; label)
    lines!(ax_θ, run.hours, run.θ_above_residual; label)
end
axislegend(ax_lhf)
save(joinpath(outdir, "soil_skin_temperature_column.png"), fig)

isnothing(with_bisection.failure) ||
    error("The column failed with the bisection: $(with_bisection.failure)")
minimum(with_bisection.θ_above_residual) > 0 || error(
    "With the bisection, the top cell dried below its residual water content",
)
if !isnothing(without_bisection) &&
   isnothing(without_bisection.failure) &&
   minimum(without_bisection.θ_above_residual) > 0
    @warn "Without the bisection, the column no longer fails"
end

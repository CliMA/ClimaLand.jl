"""
Run the calibration forward model at the prior mean parameters for a short
period, to verify the model setup before running the full calibration.

Configuration via environment variables:
    NEON_SITE_ID     — NEON site ID (default: "NEON-srer")
    NEON_SPINUP_DAYS — Number of spinup days (default: 20)

Usage:
    julia --project=.buildkite experiments/calibrate_neon/run_prior_mean.jl
"""

import ClimaLand
import ClimaLand.Parameters as LP
import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaLand.Simulations: LandSimulation, solve!
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
#using ClimaLand.Canopy.PlantHydraulics
using ClimaLand.Snow
using ClimaCore
using ClimaDiagnostics
using ClimaUtilities

import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.Utils: searchsortednearest, linear_interpolation
import Interpolations
import ClimaUtilities.TimeManager: date
import ClimaParams as CP
using Insolation

using Dates
using Statistics
using CairoMakie
using CSV
using DataFrames
CairoMakie.activate!()

include(    joinpath(pkgdir(ClimaLand), "experiments/calibrate_neon/site_metadata.jl"))

const FT = Float64
try
    global climaland_dir = pkgdir(ClimaLand)
catch err
    @warn "Using the following CLimaland dir:."
    println(pwd())
end
#const SITE_ID = get(ENV, "NEON_SITE_ID", "NEON-srer")
#const SPINUP_DAYS = parse(Int, get(ENV, "NEON_SPINUP_DAYS", "20"))
const DT = Float64(180) #(450)
const N_ITERATIONS = parse(Int, get(ENV, "NEON_N_ITERATIONS", "10"))
const Caldepthnum = get(ENV, "CALL_DEPTH", "0.00")
const outfolder = get(ENV, "outfolder", "output")

outputpath = get(ENV, "CALL_OUTPUT_PATH", "/kiwi-data/Data/groupMembers/evametz/ClimaLand_Output/Neon_siteruns/$(SITE_ID)/")

# Pipeline mode: when CALL_FIGURES_DIR is set (by start_calibration_pipeline.jl),
# write figures into that folder inside the calibration output and read the
# optimized parameters from the just-completed calibration's
# final_parameter_means.txt. When unset, behave as before (standalone use).
const PIPELINE_MODE = haskey(ENV, "CALL_FIGURES_DIR")

#replace $(N_ITERATIONS)-It with "prior_mean" for the prior mean run output path
Prior_filepath = if PIPELINE_MODE
    get(ENV, "CALL_PRIOR_FILEPATH",
        joinpath(get(ENV, "CALL_OUTPUT_1", ""), "final_parameter_means.txt"))
else
    ""#"$(outputpath)/$(outfolder)/final_parameter_means.txt"
end

if !PIPELINE_MODE
    outputpath = replace(outputpath, "$(N_ITERATIONS)-It" => "prior_mean")
end
#-----HELPER FUNCTIONS

function find_available_dir(base_dir)
    """Find an available directory by appending _1, _2, etc. if base_dir exists."""
    if !isdir(base_dir)
        return base_dir
    end
    i = 1
    while isdir("$(base_dir)_$i")
        i += 1
    end
    return "$(base_dir)_$i"
end

#------

time_offset = 0
(site_start_date, site_stop_date) =
    FluxnetSimulations.get_data_dates(SITE_ID, time_offset)
#overwrite
start_date =
    DateTime(get(ENV, "NEON_START_DATE", string(Date(site_start_date))))
stop_date = DateTime(get(ENV, "NEON_STOP_DATE", string(Date(site_stop_date))))

spinup_date = start_date + Day(SPINUP_DAYS)

if PIPELINE_MODE
    OUTPUT_DIR = ENV["CALL_FIGURES_DIR"]
else
    output_base = joinpath(outputpath, "prior_runs/prior_mean")
    OUTPUT_DIR = find_available_dir(output_base)
end
mkpath(OUTPUT_DIR)
outpath = joinpath(OUTPUT_DIR, "prior_mean_$(SITE_ID).png")

# copy folder with model scripts to output dir for record-keeping
scripts_src = joinpath(climaland_dir, "experiments/calibrate_neon")
scripts_dst = joinpath(OUTPUT_DIR, "model_scripts")
cp(scripts_src, scripts_dst; force=true)
scripts_src = joinpath(climaland_dir, "src/standalone/Soil/Biogeochemistry")
cp(scripts_src, joinpath(scripts_dst, "Biogeochemistry"); force=true)

time_offset = 0
metadata = _get_neon_site_metadata(SITE_ID)
lat = FT(metadata.lat)
long = FT(metadata.long)
atmos_h = FT(metadata.atmos_h)

# ── Prior mean values ─────────────────────────────────────────────────────────
const PRIOR_TOML = joinpath(OUTPUT_DIR, "prior_mean_parameters.toml")
function set_prior()
    soilCO2_reference_rate = 8.0329e-8
    soilCO2_activation_energy = 161640.0
    michaelis_constant = 0.049396
    O2_michaelis_constant = 0.035488

    # read the prior values from final_parameter_means.txt if it exists, otherwise use the defaults above
    if isfile(Prior_filepath)
        println("Reading prior mean parameters from $Prior_filepath")
        for line in eachline(Prior_filepath)
            m = match(r"^\s*(\w+)\s*=\s*([-\deE.+]+)", line)
            m === nothing && continue
            key, val = m.captures[1], parse(Float64, m.captures[2])
            if key == "soilCO2_reference_rate"
                soilCO2_reference_rate = val
            elseif key == "soilCO2_activation_energy"
                soilCO2_activation_energy = val
            elseif key == "michaelis_constant"
                michaelis_constant = val
            elseif key == "O2_michaelis_constant"
                O2_michaelis_constant = val
            end
        end
        println("  soilCO2_reference_rate    = $soilCO2_reference_rate")
        println("  soilCO2_activation_energy = $soilCO2_activation_energy")
        println("  michaelis_constant        = $michaelis_constant")
        println("  O2_michaelis_constant     = $O2_michaelis_constant")
    else
        println("Prior mean parameters file not found at $Prior_filepath, using default values.")
    end
    return soilCO2_reference_rate, soilCO2_activation_energy, michaelis_constant, O2_michaelis_constant
end
soilCO2_reference_rate, soilCO2_activation_energy, michaelis_constant, O2_michaelis_constant = set_prior()
println("  soilCO2_reference_rate    = $soilCO2_reference_rate")
println("  soilCO2_activation_energy = $soilCO2_activation_energy")
println("  michaelis_constant        = $michaelis_constant")
println("  O2_michaelis_constant     = $O2_michaelis_constant")

open(PRIOR_TOML, "w") do io
    write(io, """
[soilCO2_reference_rate]
value = $(soilCO2_reference_rate)
type = "float"
used_in = ["Land"]

[soilCO2_activation_energy]
value = $(soilCO2_activation_energy)
type = "float"
used_in = ["Land"]

[michaelis_constant]
value = $(michaelis_constant)
type = "float"
used_in = ["Land"]

[O2_michaelis_constant]
value = $(O2_michaelis_constant)
type = "float"
used_in = ["Land"]
""")
end
println("Wrote prior mean TOML: $PRIOR_TOML")

# ── Setup ─────────────────────────────────────────────────────────────────────
site_ID_val = FluxnetSimulations.replace_hyphen(SITE_ID)

#(; time_offset, lat, long) =
#    FluxnetSimulations.get_location(FT, Val(site_ID_val))
#(start_date, stop_date) =
#    FluxnetSimulations.get_data_dates(SITE_ID, time_offset)
#spinup_date = start_date + Day(SPINUP_DAYS)

println("Site: $SITE_ID")
println("Simulating $start_date → $stop_date (spinup until $spinup_date)")

# Domain
dz_bottom = FT(2) #FT(1.5),
dz_top = FT(0.038)
dz_tuple = (dz_bottom, dz_top)
nelements = 10#24
zmin = FT(-6.2)
zmax = FT(0)
# Domain
#(; dz_tuple, nelements, zmin, zmax) =
#    FluxnetSimulations.get_domain_info(FT, Val(site_ID_val))
#(; atmos_h) =
#    FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val))

land_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
    dz_tuple = dz_tuple,
    longlat = (long, lat),
)
canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)
surface_space = land_domain.space.surface

# Determine target layer for soil CO₂ extraction
z_field = ClimaCore.Fields.coordinate_field(land_domain.space.subsurface).z
z_vals = parent(z_field)[:, 1]
target_depth = FT(-1 * parse(Float64, Caldepthnum))
target_layer = argmin(abs.(z_vals .- target_depth))
println("Target layer: $target_layer (z = $(z_vals[target_layer]) m)")

# Base and calibrated TOML
toml_dict_base = LP.create_toml_dict(FT)
toml_dict = LP.create_toml_dict(FT; override_files = [PRIOR_TOML])
#toml_dict.data["soil_C_substrate_diffusivity"]["value"] = FT(1)  # default 3.17
#toml_dict.data["CO2_diffusion_coefficient"]["value"] = FT(3.0e-5)  # default 1.39e-5 

# ERA5 forcing
(; atmos, radiation) = FluxnetSimulations.prescribed_forcing_fluxnet(
    SITE_ID,
    lat,
    long,
    time_offset,
    atmos_h,
    start_date,
    toml_dict_base,
    FT,
)

# Custom canopy aerodynamic coefficients
toml_dict_base.data["canopy_d_coeff"]["value"] = FT(0.67)
toml_dict_base.data["canopy_z_0b_coeff"]["value"] = FT(0.013)
toml_dict_base.data["canopy_z_0m_coeff"]["value"] = FT(0.13)

# MODIS LAI
#LAI = ClimaLand.Canopy.prescribed_lai_modis(surface_space, start_date, stop_date)
LAI = ClimaLand.Canopy.prescribed_climatological_lai_modis(surface_space)


# ── Build model ──────────────────────────────────────────────────────────────
prognostic_land_components = (:canopy, :snow, :soil, :soilco2)
forcing = (; atmos, radiation)
ground = ClimaLand.PrognosticGroundConditions{FT}()
canopy_forcing = (; atmos, radiation, ground)

photosynthesis = PModel{FT}(land_domain, toml_dict_base)
conductance = PModelConductance{FT}(toml_dict_base)
soil_moisture_stress =
    ClimaLand.Canopy.PiecewiseMoistureStressModel{FT}(land_domain, toml_dict_base)

canopy = ClimaLand.Canopy.CanopyModel{FT}(
    canopy_domain,
    canopy_forcing,
    LAI,
    toml_dict_base;
    prognostic_land_components,
    photosynthesis,
    conductance,
    soil_moisture_stress,
)

α_snow = Snow.ZenithAngleAlbedoModel(toml_dict_base)
snow = Snow.SnowModel(
    FT,
    canopy_domain,
    forcing,
    toml_dict_base,
    DT;
    prognostic_land_components,
    α_snow,
)

land = LandModel{FT}(
    forcing,
    LAI,
    toml_dict,
    land_domain;
    prognostic_land_components,
    snow,
    canopy,
)

porosity_scale = FT(1)
land.soil.parameters.ν .*= porosity_scale
println("Porosity scaled by $porosity_scale, mean ν = $(mean(parent(land.soil.parameters.ν)))")


# ── Initial conditions with SOC profile ──────────────────────────────────────
base_set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
    SITE_ID,
    start_date,
    time_offset,
    land,
)
#=
function set_ic!(Y, p, t, model)
    base_set_ic!(Y, p, t, model)
    Y.soilco2.CO2 .= FT(6e-5)
    Y.soilco2.O2 .= FT(0.08)
    SOC_top = FT(18.0)
    SOC_bot = FT(0.5)
    τ_soc = FT(1.0 / log(SOC_top / SOC_bot))
    z = ClimaCore.Fields.coordinate_field(axes(Y.soilco2.SOC)).z
    @. Y.soilco2.SOC = SOC_bot + (SOC_top - SOC_bot) * exp(z / τ_soc)
end=#


ocd_path = ClimaLand.Artifacts.soil_grids_ocd_artifact_path()
SOC_from_artifact = SpaceVaryingInput(
    ocd_path,
    "ocd",
    land_domain.space.subsurface;
    regridder_type = :InterpolationsRegridder,
    regridder_kwargs = (;
        extrapolation_bc = (
            Interpolations.Periodic(),
            Interpolations.Flat(),
            Interpolations.Flat(),
        ),
        interpolation_method = Interpolations.Linear(),
    ),
)
#=
function set_ic!(Y, p, t, model)
    base_set_ic!(Y, p, t, model)
    Y.soilco2.CO2 .= FT(6e-5)
    Y.soilco2.O2 .= FT(0.08)
    #read csv file with depth and SOC values, then interpolate to model layers
    model_value = ClimaCore.Fields.zeros(land_domain.space.subsurface)
    data = CSV.read("/kiwi-data/Data/groupMembers/evametz/Neon/Neon_data/NEON_all_sites_estimatedOC_2cm_mean_extrapolated.csv", DataFrame)
    valid = .!ismissing.(data[!, "$(SITE_ID)_estimatedOC_kg_m3"])
    #z_bottom = minimum(parent(ClimaCore.Fields.coordinate_field(land_domain.space.subsurface).z))
    raw_z::Vector{Float64} = Float64.(data.depth[valid])
    raw_vals::Vector{Float64} = Float64.(data[valid, "$(SITE_ID)_estimatedOC_kg_m3"])
    sort_idx = sortperm(raw_z)
    # z_bottom is the most negative value → prepend so itp_z stays ascending
    #itp_z::Vector{Float64} = vcat(z_bottom, raw_z[sort_idx])
    #itp_values::Vector{Float64} = vcat(FT(0.5), raw_vals[sort_idx])
    zvalues = ClimaCore.Fields.coordinate_field(axes(Y.soilco2.SOC)).z
    #model_value .= map(zvalues) do z
    #    linear_interpolation(itp_z, itp_values, z)
    model_value .= map(zvalues) do z
        linear_interpolation(raw_z[sort_idx], raw_vals[sort_idx], z)
    end
    Y.soilco2.SOC .= model_value
end=#

#= old SOC-only set_ic! (uniform-column SWC from base_set_ic!)
function set_ic!(Y, p, t, model)
    base_set_ic!(Y, p, t, model)
    #Y.soilco2.CO2 .= FT(6e-5)
    #Y.soilco2.O2 .= FT(0.08)
    #read csv file with depth and SOC values, then interpolate to model layers
    model_value = ClimaCore.Fields.zeros(land_domain.space.subsurface)
    data = CSV.read("/kiwi-data/Data/groupMembers/evametz/Neon/Neon_data/NEON_all_sites_estimatedOC_2cm_mean.csv", DataFrame)
    valid = .!ismissing.(data[!, "$(SITE_ID)_estimatedOC_kg_m3"])
    raw_z::Vector{Float64} = Float64.(data.depth[valid])
    sort_idx = sortperm(raw_z)
    raw_vals::Vector{Float64} = Float64.(data[valid, "$(SITE_ID)_estimatedOC_kg_m3"])

    #get extrapolation values
    # get the first entry in raw_z, which is the most negative value (deepest depth)
    z_extrap_top = (raw_z[sort_idx])[1]
    SOC_extrap_top = (raw_vals[sort_idx])[1]
    SOC_extrap_bot = FT(0.5)
    z_extrap_bot = minimum(parent(ClimaCore.Fields.coordinate_field(land_domain.space.subsurface).z))

    zvalues = ClimaCore.Fields.coordinate_field(axes(Y.soilco2.SOC)).z

    alpha_soc = FT(log(SOC_extrap_top / SOC_extrap_bot) / (z_extrap_bot - z_extrap_top))

    # append the two parts to get the full profile
    model_value .= map(zvalues) do z
        if z > z_extrap_top
            linear_interpolation(raw_z[sort_idx], raw_vals[sort_idx], z)
        else
            SOC_extrap_top * exp(- alpha_soc * (z - z_extrap_top))
        end
    end
    Y.soilco2.SOC .= model_value

end
=#

# New set_ic! — SOC profile from NEON CSV + soil-moisture profile derived from
# NEON VSWCMean sensors (time- and plot-mean per depth code).
function set_ic!(Y, p, t, model)
    base_set_ic!(Y, p, t, model)

    # ── 1. SOC profile (NEON CSV with exponential extrapolation below data) ──
    soc_field = ClimaCore.Fields.zeros(land_domain.space.subsurface)
    soc_data = CSV.read(
        "/kiwi-data/Data/groupMembers/evametz/Neon/Neon_data/NEON_all_sites_estimatedOC_2cm_mean.csv",
        DataFrame,
    )
    valid_soc = .!ismissing.(soc_data[!, "$(SITE_ID)_estimatedOC_kg_m3"])
    raw_z::Vector{Float64} = Float64.(soc_data.depth[valid_soc])
    sort_idx_soc = sortperm(raw_z)
    raw_vals::Vector{Float64} =
        Float64.(soc_data[valid_soc, "$(SITE_ID)_estimatedOC_kg_m3"])

    z_extrap_top = (raw_z[sort_idx_soc])[1]
    SOC_extrap_top = (raw_vals[sort_idx_soc])[1]
    SOC_extrap_bot = FT(0.05)
    z_extrap_bot = minimum(parent(
        ClimaCore.Fields.coordinate_field(land_domain.space.subsurface).z,
    ))
    zvalues = ClimaCore.Fields.coordinate_field(axes(Y.soilco2.SOC)).z
    alpha_soc =
        FT(log(SOC_extrap_top / SOC_extrap_bot) / (z_extrap_bot - z_extrap_top))

    soc_field .= map(zvalues) do z
        if z > z_extrap_top
            linear_interpolation(raw_z[sort_idx_soc], raw_vals[sort_idx_soc], z)
        else
            SOC_extrap_top * exp(-alpha_soc * (z - z_extrap_top))
        end
    end
    Y.soilco2.SOC .= soc_field

    # ── 2. Soil moisture profile from NEON VSWCMean columns ─────────────────
    # Standard NEON soil sensor depths (m, negative = below surface).
    # 501 ≈ 6 cm, 502 ≈ 16 cm, 503 ≈ 26 cm, 504 ≈ 46 cm, 505 ≈ 66 cm,
    # 506 ≈ 86 cm, 507 ≈ 106 cm, 508 ≈ 166 cm
    neon_depths = FT[-0.06, -0.16, -0.26, -0.46, -0.66, -0.86, -1.06, -1.66]
    depth_codes = ["501", "502", "503", "504", "505", "506", "507", "508"]
    n_plots = 5

    csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(SITE_ID)
    swc_data = CSV.read(csv_path, DataFrame)
    swc_colnames = names(swc_data)

    # Time- and plot-mean SWC per depth code
    swc_per_depth = FT[]
    for code in depth_codes
        vals = Float64[]
        for plot_id in 1:n_plots
            colname = "VSWCMean_$(lpad(plot_id, 3, '0'))_$code"
            colname in swc_colnames || continue
            for v in swc_data[!, colname]
                (ismissing(v) || isnan(Float64(v))) && continue
                push!(vals, Float64(v))
            end
        end
        push!(swc_per_depth, isempty(vals) ? FT(NaN) : FT(mean(vals)))
    end

    valid_swc = .!isnan.(swc_per_depth)
    swc_z_valid = neon_depths[valid_swc]
    swc_vals_valid = swc_per_depth[valid_swc]
    sort_idx_swc = sortperm(swc_z_valid)
    swc_z_sorted = swc_z_valid[sort_idx_swc]
    swc_vals_sorted = swc_vals_valid[sort_idx_swc]

    println("NEON-derived SWC profile (z [m], θ_l):")
    for (zd, vd) in zip(swc_z_sorted, swc_vals_sorted)
        println("  z = $(round(zd, digits=3)), θ_l = $(round(vd, digits=4))")
    end

    z_top_data = swc_z_sorted[end]
    z_bot_data = swc_z_sorted[1]
    swc_top = swc_vals_sorted[end]
    swc_bot = swc_vals_sorted[1]

    z_soil = ClimaCore.Fields.coordinate_field(axes(Y.soil.ϑ_l)).z
    Y.soil.ϑ_l .= map(z_soil) do z
        if z > z_top_data
            FT(swc_top)
        elseif z < z_bot_data
            FT(swc_bot)
        else
            FT(linear_interpolation(swc_z_sorted, swc_vals_sorted, z))
        end
    end

    # Clip to physical range (θ_r + ε, ν - ε) to avoid retention-curve issues
    ν_field = land.soil.parameters.ν
    θ_r_field = land.soil.parameters.θ_r
    @. Y.soil.ϑ_l =
        clamp(Y.soil.ϑ_l, θ_r_field + FT(1e-4), ν_field - FT(1e-4))
end

#=
function set_ic!(Y, p, t, model)
    base_set_ic!(Y, p, t, model)
    Y.soilco2.CO2 .= FT(6e-5)
    Y.soilco2.O2 .= FT(0.08)
    Y.soilco2.SOC .= SOC_from_artifact
end=#

# ── Diagnostics ──────────────────────────────────────────────────────────────
output_writer = ClimaDiagnostics.Writers.DictWriter()
output_vars = [
    "swc",
    "tsoil",
    "si",
    "sco2",
    "soc",
    "hr",
    "so2",
    "sco2_ppm",
    "scd",
    "scms",
]
diags = ClimaLand.default_diagnostics(
    land,
    start_date;
    output_writer = output_writer,
    output_vars,
    reduction_period = :halfhourly,
)

simulation = LandSimulation(
    start_date,
    stop_date,
    DT,
    land;
    set_ic! = set_ic!,
    updateat = Second(DT),
    diagnostics = diags,
)

println("Running prior-mean model...")
@time solve!(simulation)
println("Done.")

# ── Extract diagnostics ──────────────────────────────────────────────────────
function get_diag_series(sim, name, layer)
    (times, data) = ClimaLand.Diagnostics.diagnostic_as_vectors(
        sim.diagnostics[1].output_writer,
        name;
        layer = layer,
    )
    model_dates = times isa Vector{DateTime} ? times : date.(times)
    df = DataFrame(datetime = model_dates, value = Float64.(data))
    df[!, :date] = Date.(df.datetime)
    df = filter(row -> row.date >= Date(spinup_date), df)
    sort!(df, :datetime)
    return df
end

function get_diag_layer(sim, name, layer)
    df = get_diag_series(sim, name, layer)
    daily = combine(groupby(df, :date), :value => mean => :daily_mean)
    sort!(daily, :date)
    return daily
end

sco2_daily = get_diag_layer(simulation, "sco2_ppm_30m_average", target_layer)
so2_daily = get_diag_layer(simulation, "so2_30m_average", target_layer)
swc_daily = get_diag_layer(simulation, "swc_30m_average", target_layer)
tsoil_daily = get_diag_layer(simulation, "tsoil_30m_average", target_layer)

# ── SOC depth profile plot ────────────────────────────────────────────────────
# SOC is time-invariant (dY.soilco2.SOC = 0), so the first timestep is enough
soc_vals = Float64[]
for layer in 1:nelements
    (_, data) = ClimaLand.Diagnostics.diagnostic_as_vectors(
        simulation.diagnostics[1].output_writer,
        "soc_30m_average";
        layer = layer,
    )
    push!(soc_vals, data[1])
end

fig_soc = Figure(size = (500, 700))
ax_soc = Axis(
    fig_soc[1, 1];
    xlabel = "SOC (kg C m⁻³)",
    ylabel = "Depth (m)",
    title  = "$SITE_ID — SOC Profile",
)
lines!(ax_soc, soc_vals, z_vals; color = :brown, linewidth = 2.0, label = "Model SOC")
scatter!(ax_soc, soc_vals, z_vals; color = :brown, markersize = 6)
axislegend(ax_soc; position = :rb, framevisible = false)

outpath_soc = joinpath(OUTPUT_DIR, "soc_profile_$(SITE_ID).png")
CairoMakie.save(outpath_soc, fig_soc)
println("Saved SOC profile plot: $outpath_soc")

#instability analyses
#=
# Halfhourly trigger terms around O2 instability
so2_hh = get_diag_series(simulation, "so2_30m_average", target_layer)
sod_hh = get_diag_series(simulation, "sod_30m_average", target_layer)
scd_hh = get_diag_series(simulation, "scd_30m_average", target_layer)
scms_hh = get_diag_series(simulation, "scms_30m_average", target_layer)
swc_hh = get_diag_series(simulation, "swc_30m_average", target_layer)

trigger_df = select(so2_hh, :datetime, :date)
trigger_df[!, :o2_f] = so2_hh.value
trigger_df[!, :o2_diffusivity] = sod_hh.value
trigger_df[!, :co2_diffusivity] = scd_hh.value
trigger_df[!, :microbe_source] = scms_hh.value
trigger_df[!, :swc] = swc_hh.value

nrows = nrow(trigger_df)
trigger_df[!, :dt_seconds] = fill(NaN, nrows)
if nrows >= 2
    trigger_df[2:end, :dt_seconds] =
        Float64.(Dates.value.(trigger_df.datetime[2:end] .- trigger_df.datetime[1:(end - 1)])) ./ 1000
end

trigger_df[!, :o2_dt] = fill(NaN, nrows)
if nrows >= 2
    Δo2 = trigger_df.o2_f[2:end] .- trigger_df.o2_f[1:(end - 1)]
    trigger_df[2:end, :o2_dt] = Δo2 ./ trigger_df.dt_seconds[2:end]
end

if target_layer > 1 && target_layer < length(z_vals)
    so2_up_hh = get_diag_series(simulation, "so2_30m_average", target_layer - 1)
    so2_dn_hh = get_diag_series(simulation, "so2_30m_average", target_layer + 1)
    dz_centered = abs(z_vals[target_layer - 1] - z_vals[target_layer + 1])
    trigger_df[!, :o2_grad_centered] =
        (so2_up_hh.value .- so2_dn_hh.value) ./ dz_centered
else
    trigger_df[!, :o2_grad_centered] = fill(NaN, nrows)
end

trigger_csv = joinpath(OUTPUT_DIR, "o2_trigger_terms_$(SITE_ID).csv")
CSV.write(trigger_csv, trigger_df)
println("Saved O2 trigger diagnostics: $trigger_csv")

valid_idx = findall(i -> isfinite(trigger_df.o2_dt[i]), eachindex(trigger_df.o2_dt))
if !isempty(valid_idx)
    sorted_idx = sort(valid_idx; by = i -> abs(trigger_df.o2_dt[i]), rev = true)
    top_n = min(10, length(sorted_idx))
    println("\nTop O2 tendency spikes (|dO2/dt|):")
    for i in sorted_idx[1:top_n]
        println(
            "  $(trigger_df.datetime[i]) | dO2/dt=$(trigger_df.o2_dt[i]) s^-1 | O2=$(trigger_df.o2_f[i]) | swc=$(trigger_df.swc[i]) | sod=$(trigger_df.o2_diffusivity[i]) | grad=$(trigger_df.o2_grad_centered[i])",
        )
    end
end
=#

# ── Load NEON observations ───────────────────────────────────────────────────
csv_path = ClimaLand.Artifacts.experiment_fluxnet_data_path(SITE_ID)
obs_df = CSV.read(csv_path, DataFrame)

if Caldepthnum == "0.02"
    ObsDepth = "501"
elseif Caldepthnum == "0.06"
    ObsDepth = "502"
else
    error("Calibration depth $Caldepthnums m not recognized. Please set CALL_DEPTH to a valid value (e.g., 0.02 for 2 cm, 0.06 for 6 cm).")
end

co2_cols_501 = [
    Symbol("soilCO2concentrationMean_001_$ObsDepth"),
    Symbol("soilCO2concentrationMean_002_$ObsDepth"),
    Symbol("soilCO2concentrationMean_003_$ObsDepth"),
    Symbol("soilCO2concentrationMean_004_$ObsDepth"),
    Symbol("soilCO2concentrationMean_005_$ObsDepth"),
]

function rowmean_skipinvalid(row, cols)
    vals = Float64[]
    for c in cols
        v = row[c]
        if !ismissing(v) && !isnan(Float64(v))
            push!(vals, Float64(v))
        end
    end
    return isempty(vals) ? NaN : mean(vals)
end

obs_df[!, :sco2_mean_501] =
    [rowmean_skipinvalid(row, co2_cols_501) for row in eachrow(obs_df)]
obs_df[!, :datetime] =
    DateTime.(string.(Int.(obs_df.timestamp_fmt)), dateformat"yyyymmddHHMM")
obs_df[!, :date] = Date.(obs_df.datetime)

obs_daily = combine(
    groupby(obs_df, :date),
    :sco2_mean_501 =>
        (x -> begin
            valid = filter(!isnan, x)
            length(valid) >= 24 ? mean(valid) : NaN
        end) => :daily_mean,
)
obs_daily = filter(row -> row.date >= Date(spinup_date), obs_daily)
obs_daily = filter(row -> !isnan(row.daily_mean), obs_daily)
sort!(obs_daily, :date)

# ── Plot ──────────────────────────────────────────────────────────────────────
fig = Figure(size=(1200, 1000))

ax1 = Axis(fig[1,1]; ylabel="Soil CO₂ (ppm)", title="$SITE_ID — Prior Mean vs NEON Obs (depth $(ObsDepth))")
lines!(ax1, sco2_daily.date, sco2_daily.daily_mean; color=:blue, linewidth=1.5, label="Prior mean model")
lines!(ax1, obs_daily.date, Float64.(obs_daily.daily_mean); color=:black, linewidth=1.5, label="NEON obs")
axislegend(ax1; position=:rt, framevisible=false)

ax2 = Axis(fig[2,1]; ylabel="SWC", xlabel="Date")
lines!(ax2, swc_daily.date, swc_daily.daily_mean; color=:blue, linewidth=1.5, label="SWC (model)")
axislegend(ax2; position=:rt, framevisible=false)

ax3 = Axis(fig[3,1]; xlabel="Date", ylabel="T (K)")
lines!(ax3, tsoil_daily.date, tsoil_daily.daily_mean; color=:red, linewidth=1.5, label="T (model)")
axislegend(ax3; position=:rt, framevisible=false)

if nrow(swc_daily) > 0
    swc_xlim = extrema(swc_daily.date)
    xlims!(ax1, swc_xlim...)
    xlims!(ax2, swc_xlim...)
    xlims!(ax3, swc_xlim...)
end

mkpath(dirname(outpath))
CairoMakie.save(outpath, fig)
println("Saved: $outpath")

# ── Figure: optimized model vs all individual NEON sensor plots ──────────────
# Three panels at the target depth (ObsDepth): CO₂, SWC, soil T.
# For each variable, all 5 NEON sensor plots (001–005) are plotted as daily
# means in light grey, alongside the optimized model line.
function _daily_from_col(df, col)
    valid = .!ismissing.(df[!, col]) .& .!isnan.(Float64.(coalesce.(df[!, col], NaN)))
    isempty(findall(valid)) && return DataFrame(date = Date[], daily_mean = Float64[])
    sub = DataFrame(date = df.date[valid], value = Float64.(df[valid, col]))
    daily = combine(groupby(sub, :date),
                    :value => (x -> begin
                        v = filter(!isnan, x)
                        length(v) >= 24 ? mean(v) : NaN
                    end) => :daily_mean)
    daily = filter(r -> !isnan(r.daily_mean), daily)
    daily = filter(r -> r.date >= Date(spinup_date), daily)
    sort!(daily, :date)
    return daily
end

n_plots_obs = 5
co2_plot_cols  = [Symbol("soilCO2concentrationMean_$(lpad(p,3,'0'))_$ObsDepth") for p in 1:n_plots_obs]
swc_plot_cols  = [Symbol("VSWCMean_$(lpad(p,3,'0'))_$ObsDepth")               for p in 1:n_plots_obs]
tsoil_plot_cols = [Symbol("soilTempMean_$(lpad(p,3,'0'))_$ObsDepth")          for p in 1:n_plots_obs]

obs_colnames = names(obs_df)

fig_opt = Figure(size = (1200, 1000))
ax_co2 = Axis(fig_opt[1, 1]; ylabel = "Soil CO₂ (ppm)",
              title = "$SITE_ID — Optimized model vs individual NEON sensors (depth $(ObsDepth))")
ax_swc = Axis(fig_opt[2, 1]; ylabel = "SWC (m³/m³)")
ax_t   = Axis(fig_opt[3, 1]; ylabel = "T (K)", xlabel = "Date")

obs_color = (:gray40, 0.6)
for (i, col) in enumerate(co2_plot_cols)
    String(col) in obs_colnames || continue
    d = _daily_from_col(obs_df, col)
    nrow(d) == 0 && continue
    lab = i == 1 ? "NEON sensor plots (001–$(lpad(n_plots_obs,3,'0')))" : nothing
    lines!(ax_co2, d.date, d.daily_mean; color = obs_color, linewidth = 1.0,
           label = lab)
end
lines!(ax_co2, sco2_daily.date, sco2_daily.daily_mean;
       color = :firebrick, linewidth = 1.8, label = "Optimized model")
axislegend(ax_co2; position = :rt, framevisible = false)

for (i, col) in enumerate(swc_plot_cols)
    String(col) in obs_colnames || continue
    d = _daily_from_col(obs_df, col)
    nrow(d) == 0 && continue
    lab = i == 1 ? "NEON sensor plots" : nothing
    lines!(ax_swc, d.date, d.daily_mean; color = obs_color, linewidth = 1.0,
           label = lab)
end
lines!(ax_swc, swc_daily.date, swc_daily.daily_mean;
       color = :firebrick, linewidth = 1.8, label = "Optimized model")
axislegend(ax_swc; position = :rt, framevisible = false)

# NEON soil temperature is in °C; convert to K to match model.
for (i, col) in enumerate(tsoil_plot_cols)
    String(col) in obs_colnames || continue
    d = _daily_from_col(obs_df, col)
    nrow(d) == 0 && continue
    d.daily_mean .+= 273.15
    lab = i == 1 ? "NEON sensor plots" : nothing
    lines!(ax_t, d.date, d.daily_mean; color = obs_color, linewidth = 1.0,
           label = lab)
end
lines!(ax_t, tsoil_daily.date, tsoil_daily.daily_mean;
       color = :firebrick, linewidth = 1.8, label = "Optimized model")
axislegend(ax_t; position = :rt, framevisible = false)

if nrow(swc_daily) > 0
    xl = extrema(swc_daily.date)
    xlims!(ax_co2, xl...)
    xlims!(ax_swc, xl...)
    xlims!(ax_t,   xl...)
end

outpath_opt = joinpath(OUTPUT_DIR, "optimized_vs_sensors_$(SITE_ID).png")
CairoMakie.save(outpath_opt, fig_opt)
println("Saved: $outpath_opt")

# Print summary stats
println("\nSoil CO₂ stats: min=$(round(minimum(sco2_daily.daily_mean), digits=1)), max=$(round(maximum(sco2_daily.daily_mean), digits=1)), mean=$(round(mean(sco2_daily.daily_mean), digits=1)) ppm")
println("Obs CO₂ stats:  min=$(round(minimum(obs_daily.daily_mean), digits=1)), max=$(round(maximum(obs_daily.daily_mean), digits=1)), mean=$(round(mean(obs_daily.daily_mean), digits=1)) ppm")
colors = cgrad(:viridis, nelements, categorical=true)

fig = Figure(size=(1200, 1000))
ax1 = Axis(fig[1,1]; ylabel="Soil CO₂ (ppm)", title="$SITE_ID — Prior Mean vs NEON Obs (depth 501)")
ax2 = Axis(fig[2,1]; ylabel="Soil O2")
ax3 = Axis(fig[3,1]; ylabel="Soil CO₂ Microbial Source")
ax4 = Axis(fig[4,1]; xlabel="Date", ylabel="Soil Water Content")
#plot O2, CO2 and cms profile
for layer in 1:nelements
    target_layer = layer
    sco2_daily = get_diag_layer(simulation, "sco2_ppm_30m_average",target_layer)
    so2_daily = get_diag_layer(simulation, "so2_30m_average", target_layer)
    scms_daily = get_diag_layer(simulation, "scms_30m_average", target_layer)
    swc_daily = get_diag_layer(simulation, "swc_30m_average", target_layer)

    lines!(ax1, sco2_daily.date, sco2_daily.daily_mean; linewidth=1.5, color=colors[target_layer], label="layer" * string(target_layer))
    lines!(ax2, so2_daily.date, so2_daily.daily_mean; linewidth=1.5, color=colors[target_layer], label="layer" * string(target_layer))
    lines!(ax3, scms_daily.date, scms_daily.daily_mean; linewidth=1.5, color=colors[target_layer], label="layer" * string(target_layer))
    lines!(ax4, swc_daily.date, swc_daily.daily_mean; linewidth=1.5, color=colors[target_layer], label="layer" * string(target_layer))

end
#axislegend(ax1; position=:rt, framevisible=false)
axislegend(ax2; position=:rt, framevisible=false)
#axislegend(ax3; position=:rt, framevisible=false)

if nrow(so2_daily) > 0
    swc_xlim = extrema(swc_daily.date)
    xlims!(ax1, swc_xlim...)
    xlims!(ax2, swc_xlim...)
    xlims!(ax3, swc_xlim...)
    xlims!(ax4, swc_xlim...)
end
outpath = joinpath(OUTPUT_DIR, "prior_O2_co2_cms_$(SITE_ID).png")

mkpath(dirname(outpath))
CairoMakie.save(outpath, fig)
println("Saved: $outpath")
#plot current SOC profile from data


#old soc depth plot
#=
#print SOC_from_artifact in a new output file for debugging
ocd_output_path = joinpath(OUTPUT_DIR, "ocd_profile.csv")
z_output_path = joinpath(OUTPUT_DIR, "z_profile.csv")
OLDocd_output_path = joinpath(OUTPUT_DIR, "ocdExp_profile.csv")

z_field = ClimaCore.Fields.coordinate_field(land_domain.space.subsurface).z
z_vals = parent(z_field)[:, 1]
#ocd_vals = evaluate!(SOC_from_artifact, z_vals)
writedlm(z_output_path, z_vals , ',')
writedlm(ocd_output_path, parent(SOC_from_artifact), ',')
SOC_top = FT(18.0)
SOC_bot = FT(0.5)
τ_soc = FT(1.0 / log(SOC_top / SOC_bot))
oldsoc = SOC_bot .+ (SOC_top - SOC_bot) .* exp.(ClimaCore.Fields.coordinate_field(land_domain.space.subsurface).z ./ τ_soc)
writedlm(OLDocd_output_path, parent(oldsoc), ',')
println("Saved SOC profile from artifact: $OLDocd_output_path")

fig = Figure(size=(600, 600))

ax1 = Axis(fig[1,1]; ylabel="Depth",xlabel="SOC (kg/m²)", title="$SITE_ID — SOC Profile Comparison")
lines!(ax1, vec(parent(SOC_from_artifact)), z_vals; color=:blue, linewidth=1.5, label="new SOC")
lines!(ax1, vec(parent(oldsoc)), z_vals; color=:black, linewidth=1.5, label="old exp SOC")
axislegend(ax1; position=:rt, framevisible=false)

outpath = joinpath(OUTPUT_DIR, "OCD_depth_$(SITE_ID).png")
mkpath(dirname(outpath))
CairoMakie.save(outpath, fig)
println("Saved: $outpath")
=#


#soil layer diagnostics for understanding O2 instability
#=
for layer in 2:nelements
    target_layer = layer
    so2_hh = get_diag_series(simulation, "so2_30m_average", target_layer)
    sod_hh = get_diag_series(simulation, "sod_30m_average", target_layer)
    scd_hh = get_diag_series(simulation, "scd_30m_average", target_layer)
    scms_hh = get_diag_series(simulation, "scms_30m_average", target_layer)
    swc_hh = get_diag_series(simulation, "swc_30m_average", target_layer)
    so2_dn_hh = get_diag_series(simulation, "so2_30m_average", target_layer - 1)
    dz = abs(z_vals[target_layer] - z_vals[target_layer - 1])
    dt_so2_hh = (so2_hh.value .- so2_dn_hh.value) ./ dz 
    dt_sco2_hh = ((get_diag_series(simulation, "sco2_ppm_30m_average", target_layer) .- get_diag_series(simulation, "sco2_ppm_30m_average", target_layer - 1))) ./ dz


    fig = Figure(size=(1200, 1000))

    ax1 = Axis(fig[1,1]; ylabel="Soil CO₂ (ppm)", title="$SITE_ID — Prior Mean vs NEON Obs (depth 501)")
    lines!(ax1, sco2_daily.date, sco2_daily.daily_mean; color=:blue, linewidth=1.5, label="Prior mean model")
    lines!(ax1, obs_daily.date, Float64.(obs_daily.daily_mean); color=:black, linewidth=1.5, label="NEON obs")
    axislegend(ax1; position=:rt, framevisible=false)

    ax2 = Axis(fig[2,1]; ylabel="SWC", xlabel="Date")
    lines!(ax2, swc_daily.date, swc_daily.daily_mean; color=:blue, linewidth=1.5, label="SWC (model)")
    axislegend(ax2; position=:rt, framevisible=false)

    ax3 = Axis(fig[3,1]; ylabel="O2 (K)")
    lines!(ax3, so2_daily.date, so2_daily.daily_mean; color=:red, linewidth=1.5, label="O2 (model)")
    axislegend(ax3; position=:rt, framevisible=false)

    ax4 = Axis(fig[4,1]; ylabel="scd")
    lines!(ax4, scd_hh.date, scd_hh.value; color=:blue, linewidth=1.5, label="Prior mean model")
    axislegend(ax4; position=:rt, framevisible=false)

    ax5 = Axis(fig[5,1]; ylabel="dt_so2_hh")
    lines!(ax5, so2_hh.date, dt_so2_hh; color=:blue, linewidth=1.5, label="dt_so2_hh")
    axislegend(ax5; position=:rt, framevisible=false)

    ax6 = Axis(fig[6,1]; ylabel="scms_hh")
    lines!(ax6, scms_hh.date, scms_hh.value; color=:red, linewidth=1.5, label="O2 (model)")
    axislegend(ax6; position=:rt, framevisible=false)

    ax7 = Axis(fig[7,1]; ylabel="sod")
    lines!(ax7, sod_hh.date, sod_hh.value; color=:blue, linewidth=1.5, label="Prior mean model")
    axislegend(ax7; position=:rt, framevisible=false)

    ax8 = Axis(fig[8,1];xlabel="Date", ylabel="dt_sco2_hh")
    lines!(ax8, sod_hh.date, dt_sco2_hh.value; color=:blue, linewidth=1.5, label="Prior mean model")
    axislegend(ax8; position=:rt, framevisible=false)



    if nrow(swc_daily) > 0
        swc_xlim = extrema(swc_daily.date)
        xlims!(ax1, swc_xlim...)
        xlims!(ax2, swc_xlim...)
        xlims!(ax3, swc_xlim...)
        xlims!(ax4, swc_xlim...)
        xlims!(ax5, swc_xlim...)
        xlims!(ax6, swc_xlim...)
        xlims!(ax7, swc_xlim...)
        xlims!(ax8, swc_xlim...)
    end

    outpath = joinpath(OUTPUT_DIR, "Soil_diagnostics_$(SITE_ID)_$(target_layer).png")
    mkpath(dirname(outpath))
    CairoMakie.save(outpath, fig)
    println("Saved: $outpath")
end
=#

# ── CO2 budget: production, emission, transport from below ───────────────────
# Extract half-hourly time series for the three budget terms

# 1. Microbial CO2 production in target layer (kg C m⁻³ s⁻¹)
scms_hh = get_diag_series(simulation, "scms_30m_average", target_layer)

# 2. CO2 emission to atmosphere: surface flux hr (mol CO2 m⁻² s⁻¹), no layer
(hr_times, hr_data) = ClimaLand.Diagnostics.diagnostic_as_vectors(
    simulation.diagnostics[1].output_writer,
    "hr_30m_average",
)
hr_dates = hr_times isa Vector{DateTime} ? hr_times : date.(hr_times)
hr_df = DataFrame(datetime = hr_dates, value = Float64.(hr_data))
hr_df[!, :date] = Date.(hr_df.datetime)
hr_df = filter(row -> row.date >= Date(spinup_date), hr_df)
sort!(hr_df, :datetime)

# 3. Diffusive CO2 transport from below: F = -D * Δ(CO2_air_eq)/Δz
#    as in the model: CO2_air_eq = CO2 / θ_eff, θ_eff ≈ ν - swc
#    Positive flux = upward transport into target_layer from below
if target_layer > 1
    scd_hh_tl  = get_diag_series(simulation, "scd_30m_average", target_layer)
    sco2_hh_tl = get_diag_series(simulation, "sco2_30m_average", target_layer)
    sco2_hh_bl = get_diag_series(simulation, "sco2_30m_average", target_layer - 1)
    swc_hh_tl  = get_diag_series(simulation, "swc_30m_average", target_layer)
    swc_hh_bl  = get_diag_series(simulation, "swc_30m_average", target_layer - 1)

    ν_tl = parent(land.soil.parameters.ν)[target_layer]
    ν_bl = parent(land.soil.parameters.ν)[target_layer - 1]

    # CO2_air_eq = CO2 / max(ν - swc, ε)  (approximation without Henry correction)
    eps_val = 1e-10
    co2_air_eq_tl = sco2_hh_tl.value ./ max.(ν_tl .- swc_hh_tl.value, eps_val)
    co2_air_eq_bl = sco2_hh_bl.value ./ max.(ν_bl .- swc_hh_bl.value, eps_val)

    dz = abs(z_vals[target_layer] - z_vals[target_layer - 1])
    transport_vals = scd_hh_tl.value .* (co2_air_eq_bl .- co2_air_eq_tl) ./ dz
    transport_df = DataFrame(
        datetime = scd_hh_tl.datetime,
        date     = scd_hh_tl.date,
        value    = transport_vals,
    )
else
    transport_df = DataFrame(
        datetime = scms_hh.datetime,
        date     = scms_hh.date,
        value    = zeros(nrow(scms_hh)),
    )
end

# ── Plot ──────────────────────────────────────────────────────────────────────
fig_budget = Figure(size = (1200, 800))

ax_prod = Axis(
    fig_budget[1, 1];
    ylabel = "CO₂ Production\n(kg C m⁻³ s⁻¹)",
    title  = "$SITE_ID — Layer CO₂ Budget (z = $(round(z_vals[target_layer], digits=3)) m)",
)
lines!(ax_prod, scms_hh.datetime, scms_hh.value; color = :darkgreen, linewidth = 1.5, label = "Microbial production")
axislegend(ax_prod; position = :rt, framevisible = false)

ax_emis = Axis(
    fig_budget[2, 1];
    ylabel = "CO₂ Emission\n(mol CO₂ m⁻² s⁻¹)",
)
lines!(ax_emis, hr_df.datetime, hr_df.value; color = :firebrick, linewidth = 1.5, label = "Emission to atmosphere (HR)")
axislegend(ax_emis; position = :rt, framevisible = false)

ax_trans = Axis(
    fig_budget[3, 1];
    ylabel = "CO₂ Transport from below\n(kg CO₂ m⁻³ s⁻¹)",
    xlabel = "Date",
)
lines!(ax_trans, transport_df.datetime, transport_df.value; color = :royalblue, linewidth = 1.5, label = "Diffusive transport from below")
hlines!(ax_trans, [0.0]; color = :gray, linewidth = 0.8, linestyle = :dash)
axislegend(ax_trans; position = :rt, framevisible = false)

if nrow(scms_hh) > 0
    t_xlim = extrema(scms_hh.datetime)
    xlims!(ax_prod,  t_xlim...)
    xlims!(ax_emis,  t_xlim...)
    xlims!(ax_trans, t_xlim...)
end

outpath_budget = joinpath(OUTPUT_DIR, "co2_budget_$(SITE_ID).png")
CairoMakie.save(outpath_budget, fig_budget)
println("Saved CO₂ budget plot: $outpath_budget")
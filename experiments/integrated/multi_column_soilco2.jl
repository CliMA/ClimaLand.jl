# # Multi-site single-column soil biogeochemistry test
#
# Runs a single column at 10 locations chosen from NaN regions of the
# global heterotrophic-respiration field, to verify the updated respiration
# parameters (lowered Ea_sx and α_sx, raised CO2 safety threshold).
#
# For each site, runs one month with ERA5 forcing and plots timeseries of
# respiration and soil biogeochemistry variables at multiple depths.

import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
import ClimaParams as CP
using Dates

using ClimaDiagnostics
using ClimaUtilities
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar

using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP
using ClimaLand.Simulations:
    LandSimulation, solve!, make_set_initial_state_from_file
using CairoMakie

const FT = Float64

# --- 10 locations in NaN regions of hr_globe.png ---
# (name, lon, lat, approximate-maxLAI)
sites = [
    ("boreal_canada", FT(-110.0), FT(55.0), FT(3.0)),
    ("alaska_interior", FT(-150.0), FT(65.0), FT(2.5)),
    ("western_siberia", FT(80.0), FT(60.0), FT(3.0)),
    ("central_europe", FT(10.0), FT(51.0), FT(4.0)),
    ("sahel", FT(8.0), FT(14.0), FT(2.0)),
    ("southern_africa", FT(25.0), FT(-23.0), FT(2.5)),
    ("australian_outback", FT(135.0), FT(-25.0), FT(2.0)),
    ("central_india", FT(78.0), FT(23.0), FT(3.5)),
    ("northern_china", FT(115.0), FT(40.0), FT(3.0)),
    ("central_us", FT(-98.0), FT(39.0), FT(3.5)),
]

# --- Time setup (same for all sites) ---
start_date = DateTime("2008-03-01")
stop_date = DateTime("2010-03-01")  # 2 years (year 1 = spinup, year 2 = plotted)
spinup_days = 365.0
Δt = 450.0  # 7.5 minutes

toml_dict = LP.create_toml_dict(FT)
context = ClimaComms.context()

# Domain vertical setup matches global run
zmin = FT(-15.0)
zmax = FT(0.0)
nelements = 15
dz_tuple = (FT(3.0), FT(0.05))

savedir = joinpath(
    pkgdir(ClimaLand),
    "experiments/integrated/multi_column_soilco2_out",
)
mkpath(savedir)

# Unit conversion: mol CO2 m-2 s-1 → gC m-2 d-1
const mol_co2_to_gC_per_day = 12.011 * 86400.0

function get_ts(writer, varname; layer = nothing)
    key = "$(varname)_1d_average"
    return ClimaLand.Diagnostics.diagnostic_as_vectors(writer, key; layer)
end

function time_to_days(times)
    t0 = times[1]
    return [Float64((t - t0).counter) / 86400.0 for t in times]
end

function run_site(name, lon, lat, maxLAI)
    @info "=============================================="
    @info "Running site: $name at (lon=$lon, lat=$lat)"
    @info "=============================================="

    domain =
        Column(; zlim = (zmin, zmax), nelements, dz_tuple, longlat = (lon, lat))
    surface_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    surface_space = domain.space.surface

    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        start_date,
        stop_date,
        surface_space,
        toml_dict,
        FT;
        use_lowres_forcing = true,
        context,
    )
    forcing = (; atmos, radiation)

    LAI = ClimaLand.Canopy.prescribed_lai_modis(
        surface_space,
        start_date,
        stop_date,
    )

    prognostic_land_components = (:canopy, :snow, :soil, :soilco2)

    photosynthesis = PModel{FT}(surface_domain, toml_dict)
    conductance = PModelConductance{FT}(toml_dict)
    soil_moisture_stress =
        ClimaLand.Canopy.PiecewiseMoistureStressModel{FT}(domain, toml_dict)
    biomass = ClimaLand.Canopy.PrescribedBiomassModel{FT}(
        domain,
        LAI,
        maxLAI,
        toml_dict,
    )
    canopy = ClimaLand.Canopy.CanopyModel{FT}(
        surface_domain,
        (;
            atmos = forcing.atmos,
            radiation = forcing.radiation,
            ground = ClimaLand.PrognosticGroundConditions{FT}(),
        ),
        LAI,
        toml_dict;
        prognostic_land_components,
        photosynthesis,
        conductance,
        soil_moisture_stress,
        biomass,
    )
    snow = Snow.SnowModel(
        FT,
        surface_domain,
        forcing,
        toml_dict,
        Δt;
        prognostic_land_components,
    )
    land = LandModel{FT}(
        forcing,
        LAI,
        toml_dict,
        domain,
        Δt;
        prognostic_land_components,
        canopy,
        snow,
    )

    output_vars =
        ["ra", "hr", "soc", "sco2_ppm", "so2", "scms", "swc", "tsoil"]
    diags = ClimaLand.default_diagnostics(
        land,
        start_date;
        output_writer = ClimaDiagnostics.Writers.DictWriter(),
        output_vars,
        reduction_period = :daily,
    )
    # Default file-based IC, then override O2_f with an exponential profile:
    # 0.21 at z=0, 0.01 at z=-1 m, continuing through the bottom of the domain.
    default_set_ic! = make_set_initial_state_from_file(
        ClimaLand.Artifacts.saturated_land_ic_path(; context),
        land,
    )
    k_o2 = FT(log(21))  # so 0.21 * exp(-k*1) = 0.01
    function set_ic_with_o2_profile!(Y, p, t0, model)
        default_set_ic!(Y, p, t0, model)
        z = model.soil.domain.fields.z
        @. Y.soilco2.O2_f = FT(0.21) * exp(k_o2 * z)
    end
    simulation = LandSimulation(
        start_date,
        stop_date,
        Δt,
        land;
        diagnostics = diags,
        set_ic! = set_ic_with_o2_profile!,
    )

    try
        @time solve!(simulation)
    catch err
        @error "Site $name failed during solve!" exception =
            (err, catch_backtrace())
        return false
    end

    # --- Quick final-state summary ---
    Y = simulation._integrator.u
    p = simulation._integrator.p
    z_coords = parent(domain.fields.z)[:]
    @info "Final state at $name (top → bottom):"
    @info "  z:           $(round.(z_coords[end:-1:1]; digits=2))"
    @info "  CO2:         $(round.(parent(Y.soilco2.CO2)[end:-1:1]; sigdigits=4))"
    @info "  O2_f:        $(round.(parent(Y.soilco2.O2_f)[end:-1:1]; sigdigits=4))"
    @info "  SOC:         $(round.(parent(Y.soilco2.SOC)[end:-1:1]; sigdigits=4))"
    @info "  θ_a:         $(round.(parent(p.soilco2.θ_a)[end:-1:1]; sigdigits=4))"
    @info "  θ_eff:       $(round.(parent(p.soilco2.θ_eff)[end:-1:1]; sigdigits=4))"
    @info "  CO2_air_eq:  $(round.(parent(p.soilco2.CO2_air_eq)[end:-1:1]; sigdigits=4))"
    @info "  Sm (src):    $(round.(parent(p.soilco2.Sm)[end:-1:1]; sigdigits=4))"
    @info "  top_bc:      $(round.(parent(p.soilco2.top_bc); sigdigits=4))"
    @info "  θ_l:         $(round.(parent(Y.soil.ϑ_l)[end:-1:1]; sigdigits=4))"

    # --- Plotting (year 2 only) ---
    writer = first(diags).output_writer
    z = parent(domain.fields.z)[:]
    nlayers = length(z)
    # Three depths: top / near-surface / shallow root zone
    selected_layers = [nlayers, nlayers - 2, nlayers - 5]
    layer_labels = ["z = $(round(z[l]; digits=2)) m" for l in selected_layers]
    ndepth = length(selected_layers)
    fracs = [0.0, 0.4, 0.7]  # shade fractions: shallow → deep
    shade(base, f) = RGBf(
        (1 - f) * base[1] + f,
        (1 - f) * base[2] + f,
        (1 - f) * base[3] + f,
    )
    black_base = (0.0, 0.0, 0.0)
    green_base = (0.0, 0.75, 0.15)
    green_rgb = RGBf(green_base...)
    left_colors = [shade(black_base, fracs[i]) for i in 1:ndepth]
    right_colors = [shade(green_base, fracs[i]) for i in 1:ndepth]

    # Helper: pull year-2 slice (rebased to 0-365 days).
    function year2(times, vals)
        days = time_to_days(times)
        keep = days .>= spinup_days
        return days[keep] .- spinup_days, vals[keep]
    end

    fs = 18
    lw = 2.0
    fig = Figure(size = (1800, 1100), fontsize = fs)
    Label(
        fig[0, :],
        "$name  (lon=$lon, lat=$lat)  —  year 2";
        fontsize = fs + 6,
    )

    # Panel 1: Ra + HR
    ax1 = Axis(fig[1, 1]; ylabel = "flux (gC m⁻² d⁻¹)")
    t_ra, v_ra = year2(get_ts(writer, "ra")...)
    t_hr, v_hr = year2(get_ts(writer, "hr")...)
    lines!(
        ax1,
        t_ra,
        v_ra .* mol_co2_to_gC_per_day;
        color = :black,
        linewidth = lw,
        label = "Ra",
    )
    lines!(
        ax1,
        t_hr,
        v_hr .* mol_co2_to_gC_per_day;
        color = RGBf(0.45, 0.45, 0.45),
        linewidth = lw,
        label = "HR",
    )
    axislegend(ax1; position = :rt)

    # Panel 2: soil CO2 (left, black) + O2_f (right, green)
    ax2 = Axis(fig[2, 1]; ylabel = "soil CO₂ (ppm)")
    ax2r = Axis(
        fig[2, 1];
        ylabel = "O₂ fraction",
        yaxisposition = :right,
        ylabelcolor = green_rgb,
        yticklabelcolor = green_rgb,
        rightspinecolor = green_rgb,
    )
    hidespines!(ax2r, :l, :t, :b)
    hidexdecorations!(ax2r)
    linkxaxes!(ax2, ax2r)
    for (i, l) in enumerate(selected_layers)
        t_c, v_c = year2(get_ts(writer, "sco2_ppm"; layer = l)...)
        t_o, v_o = year2(get_ts(writer, "so2"; layer = l)...)
        lines!(
            ax2,
            t_c,
            v_c;
            color = left_colors[i],
            linewidth = lw,
            label = layer_labels[i],
        )
        lines!(ax2r, t_o, v_o; color = right_colors[i], linewidth = lw)
    end
    axislegend(ax2; position = :rt)

    # Panel 3: soil T (left, black) + soil moisture (right, green)
    ax3 = Axis(fig[3, 1]; ylabel = "soil T (K)")
    ax3r = Axis(
        fig[3, 1];
        ylabel = "θ (m³ m⁻³)",
        yaxisposition = :right,
        ylabelcolor = green_rgb,
        yticklabelcolor = green_rgb,
        rightspinecolor = green_rgb,
    )
    hidespines!(ax3r, :l, :t, :b)
    hidexdecorations!(ax3r)
    linkxaxes!(ax3, ax3r)
    for (i, l) in enumerate(selected_layers)
        t_T, v_T = year2(get_ts(writer, "tsoil"; layer = l)...)
        t_θ, v_θ = year2(get_ts(writer, "swc"; layer = l)...)
        lines!(ax3, t_T, v_T; color = left_colors[i], linewidth = lw)
        lines!(ax3r, t_θ, v_θ; color = right_colors[i], linewidth = lw)
    end

    # Panel 4: Sm (left, black) + SOC (right, green)
    ax4 = Axis(
        fig[4, 1];
        xlabel = "day of year 2",
        ylabel = "Sm (mgC m⁻³ s⁻¹)",
    )
    ax4r = Axis(
        fig[4, 1];
        ylabel = "SOC (kgC m⁻³)",
        yaxisposition = :right,
        ylabelcolor = green_rgb,
        yticklabelcolor = green_rgb,
        rightspinecolor = green_rgb,
    )
    hidespines!(ax4r, :l, :t, :b)
    hidexdecorations!(ax4r)
    linkxaxes!(ax4, ax4r)
    for (i, l) in enumerate(selected_layers)
        t_s, v_s = year2(get_ts(writer, "scms"; layer = l)...)
        t_c, v_c = year2(get_ts(writer, "soc"; layer = l)...)
        lines!(ax4, t_s, v_s .* 1e6; color = left_colors[i], linewidth = lw)
        lines!(ax4r, t_c, v_c; color = right_colors[i], linewidth = lw)
    end

    for ax in (ax1, ax2, ax3)
        hidexdecorations!(ax; grid = false)
    end
    linkxaxes!(ax1, ax2, ax3, ax4)
    xlims!(ax1, 0, 365)

    rowgap!(fig.layout, 6)
    colsize!(fig.layout, 1, Relative(1.0))
    resize_to_layout!(fig)

    fname = joinpath(savedir, "$(name)_soilco2_timeseries.png")
    save(fname, fig)
    @info "Saved plot to $fname"
    return true
end

results = Dict{String, Bool}()
for (name, lon, lat, maxLAI) in sites
    results[name] = run_site(name, lon, lat, maxLAI)
end

@info "=============================================="
@info "Summary of site runs:"
for (name, ok) in results
    @info "  $(ok ? "OK  " : "FAIL") — $name"
end
@info "Plots in: $savedir"

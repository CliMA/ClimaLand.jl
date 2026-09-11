using Serialization, Statistics, Printf
using CairoMakie
CairoMakie.activate!(type = "png")

RESULTS = get(ENV, "RESULTS", "out/results.jls")
CONTROL = get(ENV, "CONTROL", "out_control/results.jls")   # same cases run with STORM_MM=0
FIGDIR = get(ENV, "FIGDIR", "figs")
mkpath(FIGDIR)
results = deserialize(RESULTS)
control = deserialize(CONTROL)
CASES = filter(c -> haskey(results, c), split(get(ENV, "PLOT_CASES", "sand_bare,loam_bare,clay_bare,loam_grass,loam_forest"), ","))
P_MM = 50.0

soils = Dict("sand" => (; ν = 0.43, θ_r = 0.045, α = 14.5, n = 2.68),
             "loam" => (; ν = 0.43, θ_r = 0.078, α = 3.6, n = 1.56),
             "clay" => (; ν = 0.38, θ_r = 0.068, α = 0.8, n = 1.09))
plants = Dict("grass" => (; LAI = 2.0, height = 0.5, rooting_depth = 0.1), "forest" => (; LAI = 5.0, height = 2.0, rooting_depth = 0.4))
PLANT_A = parse(Float64, get(ENV, "PLANT_A", "5e-5")); PLANT_NU = 0.2
STORM_HOURS = 12
θ_vg(s, ψ) = (m = 1 - 1 / s.n; s.θ_r + (s.ν - s.θ_r) * (1 + (s.α * abs(ψ))^s.n)^(-m))

label(c) = (p = split(c, "_"); p[2] == "bare" ? "bare $(p[1])" : "$(p[1]) + $(p[2])")
soil_of(c) = String(split(c, "_")[1]); cover_of(c) = String(split(c, "_")[2])
series(r, k) = [x[1] for x in r[k].data]
# recover cell faces from centers (centers are face midpoints; top face at z = 0)
function faces_from_centers(z)
    f = zeros(length(z) + 1); f[end] = 0.0
    for i in length(z):-1:1; f[i] = 2z[i] - f[i + 1]; end
    return f
end

# ---- per-case processed data
struct CaseData
    name::String; hours::Vector{Float64}; z::Vector{Float64}; dz::Vector{Float64}
    θ::Matrix{Float64}             # nz × nt
    runoff::Vector{Float64}; drain::Vector{Float64}; soilevap::Vector{Float64}; trans::Vector{Float64}  # mm/h
    storage::Vector{Float64}       # mm, column water (θ·dz summed), per hour
    θ0::Float64; lwp::Vector{Float64}; msf::Vector{Float64}; plant_store::Vector{Float64}; rootflux::Vector{Float64}
end
function process(name, results)
    r = results[name]; z = r["z"]; f = faces_from_centers(z); dz = diff(f)
    hours = r["swc_1h_average"].times ./ 3600
    θ = reduce(hcat, r["swc_1h_average"].data)
    veg = cover_of(name) != "bare"
    runoff = series(r, "sr_1h_average") .* 3.6e6                 # m/s → mm/h
    drain = -series(r, "sdr_1h_average") .* 3.6e6
    if veg
        et = series(r, "et_1h_average") .* 3600                  # kg/m²/s → mm/h
        trans = series(r, "trans_1h_average") .* 3600
        soilevap = et .- trans
        lwp_m = series(r, "lwp_1h_average")
        lwp = lwp_m .* 9800 ./ 1e6          # m → MPa
        msf = series(r, "msf_1h_average")
        v = plants[cover_of(name)]
        plant_store = v.LAI .* v.height .* PLANT_NU .* (1 .+ PLANT_A .* lwp_m) .* 1000   # mm of water held in the plant
        rootflux = series(r, "far_1h_average") .* 3.6e6
    else
        soilevap = series(r, "et_1h_average") .* 3.6e6           # m/s → mm/h
        trans = zeros(length(hours)); lwp = fill(NaN, length(hours)); msf = fill(NaN, length(hours)); plant_store = zeros(length(hours)); rootflux = zeros(length(hours))
    end
    storage = vec(sum(θ .* dz; dims = 1)) .* 1000
    θ0 = θ_vg(soils[soil_of(name)], -2.0)
    CaseData(name, hours, z, dz, θ, runoff, drain, soilevap, trans, storage, θ0, lwp, msf, plant_store, rootflux)
end
cases = Dict(c => process(c, results) for c in CASES)
DEPTH = -minimum(faces_from_centers(cases[CASES[1]].z))   # column depth (m)
ctrl = Dict(c => process(c, control) for c in CASES)
daily(v) = [sum(v[(24i + 1):(24i + 24)]) for i in 0:(length(v) ÷ 24 - 1)]

# ---- water budget: totals for the storm run, the no-storm control, and their difference (the storm's own fate)
totals(d) = (; R = sum(d.runoff), D = sum(d.drain), E = sum(d.soilevap), T = sum(d.trans),
    ΔS = d.storage[end] - d.θ0 * sum(d.dz) * 1000, Δplant = d.plant_store[end] - d.plant_store[1])
budget = Dict{String, NamedTuple}()
println("case          run   | runoff | drainage | soil evap | transp | ΔS soil | Δplant")
for c in CASES
    a, b = totals(cases[c]), totals(ctrl[c])
    storm = (; R = a.R - b.R, D = a.D - b.D, E = a.E - b.E, T = a.T - b.T, retained = P_MM - (a.R - b.R) - (a.D - b.D) - (a.E - b.E) - (a.T - b.T))
    budget[c] = (; storm, total = a, control = b)
    for (tag, x) in (("storm", a), ("control", b))
        @printf("%-12s %-8s %6.1f %8.1f %9.1f %8.1f %8.1f %7.1f\n", c, tag, x.R, x.D, x.E, x.T, x.ΔS, x.Δplant)
    end
    @printf("%-12s %-8s %6.1f %8.1f %9.1f %8.1f   retained %.1f of %.0f\n", "", "storm−ctl", storm.R, storm.D, storm.E, storm.T, storm.retained, P_MM)
end

# ---- colors
col = (; runoff = "#4C78A8", drain = "#72B7B2", soilevap = "#E45756", trans = "#54A24B", stored = "#B8B8B8", baseline = "#D9D9D9")
θcmap = :YlGnBu; θrange = (0.0, 0.45)

# ---- Figure 1: Hovmöller for the three bare soils
begin
    fig = Figure(size = (1200, 420), fontsize = 15)
    bare = filter(c -> cover_of(c) == "bare", CASES)
    for (i, c) in enumerate(bare)
        d = cases[c]; b = budget[c].storm
        ax = Axis(fig[1, i]; title = @sprintf("%s\nran off %.1f · evaporated %.1f · retained %.1f mm", label(c), b.R, b.E, b.retained),
            xlabel = "days since the storm", ylabel = i == 1 ? "depth (m)" : "", titlesize = 14)
        heatmap!(ax, d.hours ./ 24, faces_from_centers(d.z), permutedims(d.θ); colormap = θcmap, colorrange = θrange)
        i > 1 && hideydecorations!(ax; grid = false)
    end
    Colorbar(fig[1, length(bare) + 1]; colormap = θcmap, colorrange = θrange, label = "soil water content (m³ water / m³ soil)")
    save(joinpath(FIGDIR, "hovmoller.png"), fig; px_per_unit = 2)
end

# ---- Figure 2: where did the storm's water go? (storm run minus no-storm control)
begin
    fig = Figure(size = (900, 480), fontsize = 15)
    ax = Axis(fig[1, 1]; xlabel = "mm of water, 30 days after the storm", yticks = (1:length(CASES), label.(CASES)),
        title = "Where did the storm's 50 mm go?", yreversed = true)
    keys_ = (:R, :D, :E, :T, :retained)
    names_ = ("surface runoff", @sprintf("drained below %.0f m", DEPTH), "evaporated from soil", "transpired by plants", "still in the soil")
    colors_ = (col.runoff, col.drain, col.soilevap, col.trans, col.stored)
    for (i, c) in enumerate(CASES)
        b = budget[c].storm; x0 = 0.0
        for (k, cc) in zip(keys_, colors_)
            v = getproperty(b, k)
            v > 0.3 && barplot!(ax, [i], [v]; offset = x0, direction = :x, color = cc, strokewidth = 0.5, strokecolor = :white)
            x0 += max(v, 0.0)
        end
    end
    xlims!(ax, 0, P_MM)
    elems = [PolyElement(color = cc) for cc in colors_]
    Legend(fig[2, 1], elems, collect(names_); orientation = :horizontal, nbanks = 1, framevisible = false)
    save(joinpath(FIGDIR, "budget.png"), fig; px_per_unit = 2)
end

# ---- Figure 3: daily evapotranspiration and the plants' "straw"
begin
    fig = Figure(size = (1000, 640), fontsize = 15)
    ax1 = Axis(fig[1, 1]; ylabel = "water returned to the air (mm / day)", title = "Evaporation + transpiration after the storm", xlabel = "")
    cs = Dict("sand_bare" => "#E0A458", "loam_bare" => "#E45756", "clay_bare" => "#8C564B", "loam_grass" => "#9BD770", "loam_forest" => "#2E7D32")
    for c in CASES
        d = cases[c]; e = daily(d.soilevap .+ d.trans)
        lines!(ax1, 1:length(e), e; label = label(c), color = get(cs, c, :black), linewidth = 3)
    end
    axislegend(ax1; position = :rt, framevisible = false)
    ax2 = Axis(fig[2, 1]; ylabel = "leaf water potential (MPa)", xlabel = "days since the storm")
    ax3 = Axis(fig[2, 1]; ylabel = "moisture-stress factor (0–1)", yaxisposition = :right)
    hidespines!(ax3); hidexdecorations!(ax3)
    for c in filter(c -> cover_of(c) != "bare", CASES)
        d = cases[c]
        lwp_daymin = [minimum(d.lwp[(24i + 1):(24i + 24)]) for i in 0:(length(d.lwp) ÷ 24 - 1)]
        msf_daymean = [mean(d.msf[(24i + 1):(24i + 24)]) for i in 0:(length(d.msf) ÷ 24 - 1)]
        lines!(ax2, 1:length(lwp_daymin), lwp_daymin; color = get(cs, c, :black), linewidth = 3, label = "$(label(c)): leaf water potential, daily minimum")
        lines!(ax3, 1:length(msf_daymean), msf_daymean; color = get(cs, c, :black), linewidth = 3, linestyle = :dash, label = "$(label(c)): moisture-stress factor, daily mean")
    end
    ylims!(ax3, -0.05, 1.05)
    axislegend(ax2; position = :lb, framevisible = false, labelsize = 12)
    axislegend(ax3; position = :rt, framevisible = false, labelsize = 12)
    linkxaxes!(ax1, ax2)
    save(joinpath(FIGDIR, "et_timeseries.png"), fig; px_per_unit = 2)
end

# ---- Figure 4: the animation
# Each panel: the soil column colored by water content (storm run), and the fate of the storm's own water as
# cumulative bars (storm run minus no-storm control): ran off, evaporated, transpired above the surface, drained
# below it, and a gauge inside the column for what is still stored, all on one scale.
# Root lengths mark the 20th to 95th percentiles of the exponential root profile.
root_quantiles(rd) = [min(-rd * log(1 - q), DEPTH) for q in (0.2, 0.4, 0.6, 0.8, 0.95)]
function draw_plant!(ax, cover, rooting_depth)
    if cover == "forest"
        x = 0.8
        poly!(ax, Rect(x - 0.04, 0.0, 0.08, 0.45); color = "#6D4C41")
        for (cy, rr) in ((0.55, 0.22), (0.73, 0.19), (0.87, 0.15))
            poly!(ax, Circle(Point2f(x, cy), rr); color = "#2E7D32")
        end
        for (L, dx) in zip(root_quantiles(rooting_depth), (0.12, -0.14, 0.06, -0.05, 0.0))
            lines!(ax, [x, x + dx], [0.0, -L]; color = "#8D6E63", linewidth = 2)
        end
    elseif cover == "grass"
        for x in 0.6:0.05:0.95
            lines!(ax, [x, x + 0.03 * sign(0.78 - x)], [0.0, 0.16 + 0.05 * sin(20x)]; color = "#7CB342", linewidth = 3)
        end
        qs = root_quantiles(rooting_depth)
        for (x, L) in zip(0.62:0.04:0.94, repeat(qs, 2))
            lines!(ax, [x, x + 0.01], [0.0, -L]; color = "#8D6E63", linewidth = 1.5)
        end
    end
end

const BARS = ((:runoff, "ran off", col.runoff, 0.03), (:soilevap, "evaporated", col.soilevap, 0.22), (:trans, "transpired", col.trans, 0.41))
const BARW = 0.12
function build_animation()
    fig = Figure(size = (1250, 700), fontsize = 15, backgroundcolor = :white)
    title = Observable("")
    Label(fig[0, 1:length(CASES)], title; fontsize = 22, font = :bold, tellwidth = false)
    θobs = Dict{String, Observable{Matrix{Float64}}}()
    bars = Dict{Tuple{String, Symbol}, Observable{Vector{Point2f}}}(); blabel = Dict{Tuple{String, Symbol}, Observable{String}}()
    blabelpos = Dict{Tuple{String, Symbol}, Observable{Float64}}()
    drainbox = Dict{String, Observable{Vector{Point2f}}}(); dlabel = Dict{String, Observable{String}}()
    stored = Dict{String, Observable{Vector{Point2f}}}(); slabel = Dict{String, Observable{String}}(); spos = Dict{String, Observable{Point2f}}()
    rain = Observable(Point2f[])
    scale = 0.45 / P_MM    # 50 mm of water → 0.45 m on the drawing
    zb = -DEPTH - 0.05
    for (i, c) in enumerate(CASES)
        d = cases[c]
        ax = Axis(fig[1, i]; title = label(c), aspect = DataAspect(), limits = ((-0.05, 1.05), (-DEPTH - 0.75, 1.35)))
        hidedecorations!(ax); hidespines!(ax)
        θobs[c] = Observable(reshape(d.θ[:, 1], 1, :))
        heatmap!(ax, [0.0, 1.0], faces_from_centers(d.z), θobs[c]; colormap = θcmap, colorrange = θrange)
        lines!(ax, [0, 1, 1, 0, 0], [0, 0, -DEPTH, -DEPTH, 0]; color = :black, linewidth = 1)
        for (k, name, cc, x0) in BARS
            (k == :trans && cover_of(c) == "bare") && continue
            bars[(c, k)] = Observable(Point2f[(x0, 0), (x0 + BARW, 0), (x0 + BARW, 0), (x0, 0)])
            poly!(ax, bars[(c, k)]; color = cc, strokecolor = :white, strokewidth = 0.5)
            blabel[(c, k)] = Observable(""); blabelpos[(c, k)] = Observable(0.03)
            text!(ax, x0 + BARW / 2, blabelpos[(c, k)]; text = blabel[(c, k)], align = (:center, :bottom), fontsize = 11, color = cc)
        end
        drainbox[c] = Observable(Point2f[(0.2, zb), (0.8, zb), (0.8, zb), (0.2, zb)])
        poly!(ax, drainbox[c]; color = ("#72B7B2", 0.9))
        dlabel[c] = Observable("")
        text!(ax, 0.5, lift(p -> p[3][2] - 0.03, drainbox[c]), text = dlabel[c], align = (:center, :top), fontsize = 12, color = "#1F5F5C")
        # stored: the storm's water still in the soil, a gauge inside the column from a mid-depth baseline
        z0 = -DEPTH / 2
        lines!(ax, [0.04, 0.16], [z0, z0]; color = :white, linewidth = 1.5)
        stored[c] = Observable(Point2f[(0.04, z0), (0.16, z0), (0.16, z0), (0.04, z0)])
        poly!(ax, stored[c]; color = col.stored, strokecolor = :white, strokewidth = 1)
        spos[c] = Observable(Point2f(0.19, z0)); slabel[c] = Observable("")
        poly!(ax, lift(p -> Rect(p[1] - 0.01, p[2] - 0.045, 0.4, 0.09), spos[c]); color = (:white, 0.8))
        text!(ax, spos[c]; text = slabel[c], align = (:left, :center), fontsize = 11, color = :black)
        scatter!(ax, rain; color = "#4C78A8", marker = :vline, markersize = 12)
        draw_plant!(ax, cover_of(c), cover_of(c) == "bare" ? 0.0 : plants[cover_of(c)].rooting_depth)
    end
    Colorbar(fig[1, length(CASES) + 1]; colormap = θcmap, colorrange = θrange, label = "soil water content (m³/m³)")
    Label(fig[2, 1:length(CASES)], @sprintf("%.0f m of soil per column, colored by water content. Bars follow the storm's own 50 mm (this run minus the same column without the storm), cumulative and on one scale: ran off (blue), evaporated (red), transpired (green), drained out of the bottom (teal), still stored in the soil (gray). Roots mark the 20th–95th percentiles of the root profile.", DEPTH);
        fontsize = 11, tellwidth = false, word_wrap = true)
    hours = cases[CASES[1]].hours
    function update!(k)
        h = hours[k]; day = floor(Int, (h - 1e-9) / 24); hod = h - 24day
        title[] = @sprintf("Day %d, %02d:00 %s", day + 1, round(Int, hod), k <= STORM_HOURS ? "— raining (50 mm in 12 h)" : "")
        rain[] = k <= STORM_HOURS ? [Point2f(rand(), 0.3 + rand()) for _ in 1:60] : Point2f[]
        for c in CASES
            d = cases[c]
            θobs[c][] = reshape(d.θ[:, k], 1, :)
            hprev = -1.0
            e = ctrl[c]
            for (kk, name, cc, x0) in BARS
                haskey(bars, (c, kk)) || continue
                v = sum(getproperty(d, kk)[1:k]) - sum(getproperty(e, kk)[1:k]); h = max(v, 0.0) * scale
                lift_ = abs(h - hprev) < 0.07 ? 0.07 : 0.0   # stagger labels of neighbours at similar heights
                bars[(c, kk)][] = Point2f[(x0, 0), (x0 + BARW, 0), (x0 + BARW, h), (x0, h)]
                blabelpos[(c, kk)][] = h + lift_ + 0.03
                blabel[(c, kk)][] = v > 0.5 ? @sprintf("%.0f mm", v) : ""
                hprev = h
            end
            D = max(sum(d.drain[1:k]) - sum(e.drain[1:k]), 0.0)
            drainbox[c][] = Point2f[(0.2, zb), (0.8, zb), (0.8, zb - D * scale), (0.2, zb - D * scale)]
            dlabel[c][] = @sprintf("drained %.1f mm", D)
            ΔS = max(d.storage[k] - e.storage[k], 0.0); z0 = -DEPTH / 2
            stored[c][] = Point2f[(0.04, z0), (0.16, z0), (0.16, z0 + ΔS * scale), (0.04, z0 + ΔS * scale)]
            spos[c][] = Point2f(0.19, z0 + ΔS * scale); slabel[c][] = @sprintf("stored %+.0f mm", ΔS)
        end
    end
    return fig, update!
end

function animate(path; framerate = 12)
    fig, update! = build_animation()
    hours = cases[CASES[1]].hours
    frames = vcat(1:1:48, 54:6:length(hours))
    record(fig, path, frames; framerate) do k
        update!(k)
    end
end
if haskey(ENV, "FRAME")   # render single frames (hour indices, comma-separated) instead of the movies
    fig, update! = build_animation()
    for k in parse.(Int, split(ENV["FRAME"], ","))
        update!(k); save(joinpath(FIGDIR, "frame_$(k).png"), fig; px_per_unit = 2)
    end
elseif get(ENV, "ANIMATE", "true") == "true"
    animate(joinpath(FIGDIR, "rain_columns.gif"))
    animate(joinpath(FIGDIR, "rain_columns.mp4"))
end
println("figures written to ", FIGDIR)

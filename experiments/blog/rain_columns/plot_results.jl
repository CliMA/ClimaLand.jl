using Serialization, Statistics, Printf
using CairoMakie
CairoMakie.activate!(type = "png")

RESULTS = get(ENV, "RESULTS", "full/results.jls")
FIGDIR = get(ENV, "FIGDIR", "figs")
mkpath(FIGDIR)
results = deserialize(RESULTS)
CASES = filter(c -> haskey(results, c), split(get(ENV, "PLOT_CASES", "sand_bare,loam_bare,clay_bare,loam_grass,loam_forest"), ","))
P_MM = 50.0

soils = Dict("sand" => (; ν = 0.43, θ_r = 0.045, α = 14.5, n = 2.68),
             "loam" => (; ν = 0.43, θ_r = 0.078, α = 3.6, n = 1.56),
             "clay" => (; ν = 0.38, θ_r = 0.068, α = 0.8, n = 1.09))
plants = Dict("grass" => (; LAI = 2.0, height = 0.5, rooting_depth = 0.3), "forest" => (; LAI = 5.0, height = 2.0, rooting_depth = 1.0))
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
function process(name)
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
cases = Dict(c => process(c) for c in CASES)
daily(v) = [sum(v[(24i + 1):(24i + 24)]) for i in 0:(length(v) ÷ 24 - 1)]

# ---- water budget table
println("case | runoff | drainage | soil evap | transpiration | root uptake | Δplant | ΔS(model) | ΔS(residual)")
budget = Dict{String, NamedTuple}()
for c in CASES
    d = cases[c]
    R, D, E, T = sum(d.runoff), sum(d.drain), sum(d.soilevap), sum(d.trans)
    S0 = d.θ0 * sum(d.dz) * 1000
    ΔS_model = d.storage[end] - S0
    ΔS_resid = P_MM - R - D - E - T
    Δplant = d.plant_store[end] - d.plant_store[1]; U = sum(d.rootflux)
    budget[c] = (; R, D, E, T, ΔS = ΔS_resid, ΔS_model, Δplant)
    @printf("%-12s %6.1f %8.1f %9.1f %13.1f %11.1f %7.1f %10.1f %12.1f\n", c, R, D, E, T, U, Δplant, ΔS_model, ΔS_resid)
end

# ---- colors
col = (; runoff = "#4C78A8", drain = "#72B7B2", soilevap = "#E45756", trans = "#54A24B", stored = "#B8B8B8", deficit = "#F58518")
θcmap = :YlGnBu; θrange = (0.0, 0.45)

# ---- Figure 1: Hovmöller for the three bare soils
begin
    fig = Figure(size = (1200, 420), fontsize = 15)
    bare = filter(c -> cover_of(c) == "bare", CASES)
    for (i, c) in enumerate(bare)
        d = cases[c]; b = budget[c]
        ax = Axis(fig[1, i]; title = @sprintf("%s\nrunoff %.0f mm · drained %.0f mm · evaporated %.0f mm", label(c), b.R, b.D, b.E),
            xlabel = "days since the storm", ylabel = i == 1 ? "depth (m)" : "", titlesize = 14)
        heatmap!(ax, d.hours ./ 24, faces_from_centers(d.z), permutedims(d.θ); colormap = θcmap, colorrange = θrange)
        i > 1 && hideydecorations!(ax; grid = false)
    end
    Colorbar(fig[1, length(bare) + 1]; colormap = θcmap, colorrange = θrange, label = "soil water content (m³ water / m³ soil)")
    save(joinpath(FIGDIR, "hovmoller.png"), fig; px_per_unit = 2)
end

# ---- Figure 2: where did the 50 mm go?
begin
    fig = Figure(size = (1000, 440), fontsize = 15)
    ax = Axis(fig[1, 1]; xlabel = "mm of water (30 days after a 50 mm storm)", yticks = (1:length(CASES), label.(CASES)),
        title = "Where did the rain go?", yreversed = true)
    keys_ = (:runoff, :drain, :soilevap, :trans, :stored)
    names_ = ("surface runoff", "drained below 2 m", "evaporated from soil", "transpired by plants", "still stored in soil")
    for (i, c) in enumerate(CASES)
        b = budget[c]
        vals = [b.R, b.D, b.E, b.T, max(b.ΔS, 0.0)]
        x0 = 0.0
        for (k, v) in zip(keys_, vals)
            v > 0 && barplot!(ax, [i], [v]; offset = x0, direction = :x, color = getproperty(col, k), strokewidth = 0.5, strokecolor = :white)
            x0 += v
        end
        if b.ΔS < 0   # plants spent more than the storm delivered
            barplot!(ax, [i], [b.ΔS]; direction = :x, color = (col.deficit, 0.9), strokewidth = 0.5, strokecolor = :white)
            lbl = b.Δplant < -1 ? @sprintf("%.0f mm from pre-storm storage\n(%.0f mm of it from the plants' own tissue)", -b.ΔS, -b.Δplant) : @sprintf("%.0f mm from pre-storm\nsoil water", -b.ΔS)
            text!(ax, b.ΔS - 1.5, i; text = lbl, align = (:right, :center), fontsize = 11)
        end
    end
    vlines!(ax, [0.0]; color = :black); vlines!(ax, [P_MM]; color = :black, linestyle = :dash)
    xlims!(ax, -145, 135)
    text!(ax, P_MM + 0.5, 0.45; text = "the 50 mm storm", align = (:left, :center), fontsize = 11)
    elems = [PolyElement(color = getproperty(col, k)) for k in keys_]
    push!(elems, PolyElement(color = col.deficit)); 
    Legend(fig[2, 1], elems, [names_..., "drawn from pre-storm storage"]; orientation = :horizontal, nbanks = 2, framevisible = false)
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
    ax3 = Axis(fig[2, 1]; ylabel = "stomatal opening factor (0–1)", yaxisposition = :right)
    hidespines!(ax3); hidexdecorations!(ax3)
    for c in filter(c -> cover_of(c) != "bare", CASES)
        d = cases[c]
        lwp_daymin = [minimum(d.lwp[(24i + 1):(24i + 24)]) for i in 0:(length(d.lwp) ÷ 24 - 1)]
        msf_daymean = [mean(d.msf[(24i + 1):(24i + 24)]) for i in 0:(length(d.msf) ÷ 24 - 1)]
        lines!(ax2, 1:length(lwp_daymin), lwp_daymin; color = get(cs, c, :black), linewidth = 3, label = "$(label(c)): daily minimum leaf water potential")
        lines!(ax3, 1:length(msf_daymean), msf_daymean; color = get(cs, c, :black), linewidth = 3, linestyle = :dash, label = "$(label(c)): stomatal factor")
    end
    ylims!(ax3, -0.05, 1.05)
    axislegend(ax2; position = :lb, framevisible = false, labelsize = 12)
    axislegend(ax3; position = :rt, framevisible = false, labelsize = 12)
    linkxaxes!(ax1, ax2)
    save(joinpath(FIGDIR, "et_timeseries.png"), fig; px_per_unit = 2)
end

# ---- Figure 4: the animation
function draw_plant!(ax, cover, rooting_depth)
    if cover == "forest"
        poly!(ax, Rect(0.46, 0.0, 0.08, 0.45); color = "#6D4C41")
        for (cy, rr) in ((0.55, 0.30), (0.75, 0.26), (0.92, 0.2))
            poly!(ax, Circle(Point2f(0.5, cy), rr); color = "#2E7D32")
        end
    elseif cover == "grass"
        for x in 0.1:0.06:0.9
            lines!(ax, [x, x + 0.03 * sign(0.5 - x)], [0.0, 0.18 + 0.05 * sin(20x)]; color = "#7CB342", linewidth = 3)
        end
    end
    if cover != "bare"   # roots to scale: most roots above the rooting-depth parameter
        for (x0, x1, f) in ((0.5, 0.5, 1.0), (0.5, 0.25, 0.65), (0.5, 0.75, 0.7), (0.5, 0.35, 0.4), (0.5, 0.65, 0.5))
            lines!(ax, [x0, x1], [0.0, -rooting_depth * f * 1.5]; color = "#8D6E63", linewidth = 2)
        end
    end
end

function animate(path; framerate = 12)
    fig = Figure(size = (1250, 700), fontsize = 15, backgroundcolor = :white)
    title = Observable("")
    Label(fig[0, 1:length(CASES)], title; fontsize = 22, font = :bold, tellwidth = false)
    θobs = Dict{String, Observable{Matrix{Float64}}}()
    pond = Dict{String, Observable{Vector{Point2f}}}(); drainbox = Dict{String, Observable{Vector{Point2f}}}()
    rlabel = Dict{String, Observable{String}}(); dlabel = Dict{String, Observable{String}}(); elabel = Dict{String, Observable{String}}()
    arrow_len = Dict{String, Observable{Float64}}()
    rain = Observable(Point2f[])
    for (i, c) in enumerate(CASES)
        d = cases[c]
        ax = Axis(fig[1, i]; title = label(c), aspect = DataAspect(), limits = ((-0.05, 1.05), (-2.75, 1.35)))
        hidedecorations!(ax); hidespines!(ax)
        θobs[c] = Observable(reshape(d.θ[:, 1], 1, :))
        heatmap!(ax, [0.0, 1.0], faces_from_centers(d.z), θobs[c]; colormap = θcmap, colorrange = θrange)
        lines!(ax, [0, 1, 1, 0, 0], [0, 0, -2, -2, 0]; color = :black, linewidth = 1)
        pond[c] = Observable(Point2f[(0, 0), (1, 0), (1, 0), (0, 0)])
        poly!(ax, pond[c]; color = ("#4C78A8", 0.85))
        drainbox[c] = Observable(Point2f[(0.2, -2.05), (0.8, -2.05), (0.8, -2.05), (0.2, -2.05)])
        poly!(ax, drainbox[c]; color = ("#72B7B2", 0.9))
        rlabel[c] = Observable(""); dlabel[c] = Observable(""); elabel[c] = Observable("")
        text!(ax, 0.02, 0.02, text = rlabel[c], align = (:left, :bottom), fontsize = 12, color = "#1F3F66")
        text!(ax, 0.5, -2.08, text = dlabel[c], align = (:center, :top), fontsize = 12, color = "#1F5F5C")
        arrow_len[c] = Observable(0.0)
        lines!(ax, lift(l -> [Point2f(0.9, 0.05), Point2f(0.9, 0.05 + l)], arrow_len[c]); color = "#E45756", linewidth = 4)
        scatter!(ax, lift(l -> [Point2f(0.9, 0.05 + l)], arrow_len[c]); color = "#E45756", marker = :utriangle, markersize = 18)
        scatter!(ax, rain; color = "#4C78A8", marker = :vline, markersize = 12)
        text!(ax, 0.98, 1.3, text = elabel[c], align = (:right, :top), fontsize = 12, color = "#B03A2E")
        draw_plant!(ax, cover_of(c), cover_of(c) == "bare" ? 0.0 : plants[cover_of(c)].rooting_depth)
    end
    Colorbar(fig[1, length(CASES) + 1]; colormap = θcmap, colorrange = θrange, label = "soil water content (m³/m³)")
    Label(fig[2, 1:length(CASES)], "Each column is 2 m of soil. Blue pond = cumulative surface runoff · teal bar = cumulative drainage below 2 m · red arrow = evaporation + transpiration (24 h mean) · roots drawn to scale";
        fontsize = 12, tellwidth = false)
    hours = cases[CASES[1]].hours
    frames = vcat(1:1:48, 54:6:length(hours))
    scale = 0.45 / P_MM    # 50 mm of runoff or drainage → 0.45 m on the drawing
    record(fig, path, frames; framerate) do k
        h = hours[k]; day = floor(Int, (h - 1e-9) / 24); hod = h - 24day
        title[] = @sprintf("Day %d, %02d:00 %s", day + 1, round(Int, hod), k <= STORM_HOURS ? "— raining (50 mm in 12 h)" : "")
        rain[] = k <= STORM_HOURS ? [Point2f(rand(), 0.3 + rand()) for _ in 1:60] : Point2f[]
        for c in CASES
            d = cases[c]
            θobs[c][] = reshape(d.θ[:, k], 1, :)
            R = sum(d.runoff[1:k]) * scale; D = sum(d.drain[1:k]) * scale
            pond[c][] = Point2f[(0, 0), (1, 0), (1, R), (0, R)]
            drainbox[c][] = Point2f[(0.2, -2.05), (0.8, -2.05), (0.8, -2.05 - D), (0.2, -2.05 - D)]
            rlabel[c][] = R > 0 ? @sprintf("runoff %.0f mm", sum(d.runoff[1:k])) : ""
            dlabel[c][] = @sprintf("drained %.1f mm", sum(d.drain[1:k]))
            k0 = max(1, k - 23); et24 = sum(d.soilevap[k0:k] .+ d.trans[k0:k]) * 24 / (k - k0 + 1)
            arrow_len[c][] = 0.12 * et24
            elabel[c][] = @sprintf("%.1f mm/day\nto the air", max(et24, 0.0))
        end
    end
end
animate(joinpath(FIGDIR, "rain_columns.gif"))
animate(joinpath(FIGDIR, "rain_columns.mp4"))
println("figures written to ", FIGDIR)

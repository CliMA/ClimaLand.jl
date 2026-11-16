# The soil organic carbon (SOC) of a spun-up state against SoilGrids, on CPU: maps of the
# top metre, zonal means of the 0–1 m and 0–2 m columns, the SOC below 2 m (where
# SoilGrids has no data), and the global stock above each depth. SoilGrids is taken as
# the model reads it (`set_soilco2_SOC_from_soilgrids!` on the model's layers and land
# points), with the depth weights of `damm_choice.jl`; land north of 60°S.
#
# Usage (CPU job):
#   julia --project=.buildkite soc_spinup_figure.jl <state.hdf5> <out_dir> [<title note>]
# Writes <out_dir>/soc_spinup_1850_vs_soilgrids.png and .csv (totals by band, Pg C).
import ClimaComms
import ClimaCore
using ClimaLand
import CairoMakie as CM
using Printf

const FT = Float64
include(joinpath(@__DIR__, "carbon_budget.jl"))
include(joinpath(@__DIR__, "state_io.jl"))
const BANDS = (("global", -60, 90), ("north_50N", 50, 90), ("30N_50N", 30, 50), ("tropics_30S_30N", -30, 30), ("south_60S_30S", -60, -30))

state_path, out_dir = ARGS[1:2]
note = length(ARGS) >= 3 ? join(ARGS[3:end], " ") : basename(state_path)
mkpath(out_dir)
context = ClimaComms.SingletonCommsContext(ClimaComms.CPUSingleThreaded())
domain = ClimaLand.Domains.global_box_domain(FT; context, mask_threshold = FT(0.99))
surface = domain.space.surface
nodes, area = land_nodes_and_areas(surface)
coords = ClimaCore.Fields.coordinate_field(surface)
lat_all = vec(Array(parent(coords.lat)))
lon_all = vec(Array(parent(coords.long)))
lat, lon = lat_all[nodes], lon_all[nodes]
z = Array(parent(ClimaCore.Fields.coordinate_field(domain.space.subsurface).z))[:, 1, 1, 1, nodes]
dz = Array(parent(ClimaCore.Fields.Δz_field(domain.space.subsurface)))[:, 1, 1, 1, nodes]
nz, nl = size(z)
# thickness (m) of each layer above the depth D, SOC taken uniform within a layer
above(D) = max.(min.(-z .+ dz ./ 2, D) .- max.(-z .- dz ./ 2, 0.0), 0.0)
w1, w2, wall = above(1.0), above(2.0), dz
column(x, w) = vec(sum(x .* w; dims = 1))   # kg C m⁻²

sg_field = ClimaCore.Fields.zeros(domain.space.subsurface)
ClimaLand.Simulations.set_soilco2_SOC_from_soilgrids!(sg_field)
soilgrids = Array(parent(sg_field))[:, 1, 1, 1, nodes]
soc = read_dataset(state_path, "Y", "soilco2", "SOC")[:, 1, 1, 1, nodes]
ok = [all(isfinite, view(soc, :, i)) && all(isfinite, view(soilgrids, :, i)) for i in 1:nl]
println("state ", state_path, "\nland nodes north of 60°S ", nl, ", left out (non-finite) ", count(!, ok))

m1, m2, mall = column(soc, w1), column(soc, w2), column(soc, wall)
s1, s2 = column(soilgrids, w1), column(soilgrids, w2)
total(x, mask) = sum(x[i] * area[i] for i in 1:nl if mask[i] && ok[i]; init = 0.0) * 1e-12   # Pg C

rows = ["band,soc_0_1m_model,soc_0_1m_soilgrids,ratio_0_1m,soc_0_2m_model,soc_0_2m_soilgrids,ratio_0_2m,soc_column_model,soc_below_2m_model"]
totals = Dict{String, NTuple{5, Float64}}()
println("\nband             0-1 m model / SoilGrids (ratio) | 0-2 m model / SoilGrids (ratio) | column, below 2 m (Pg C)")
for (name, lo, hi) in BANDS
    mask = lo .<= lat .< hi
    t = (total(m1, mask), total(s1, mask), total(m2, mask), total(s2, mask), total(mall, mask))
    totals[name] = t
    @printf("%-16s %6.0f / %6.0f (%.2f) | %6.0f / %6.0f (%.2f) | %6.0f, %6.0f\n", name, t[1], t[2], t[1] / t[2], t[3], t[4], t[3] / t[4], t[5], t[5] - t[3])
    push!(rows, join((name, (@sprintf("%.2f", v) for v in (t[1], t[2], t[1] / t[2], t[3], t[4], t[3] / t[4], t[5], t[5] - t[3]))...), ","))
end
write(joinpath(out_dir, "soc_spinup_1850_vs_soilgrids.csv"), join(rows, "\n") * "\n")

# lon × lat maps from the coordinate values (1° grid)
ulon = sort(unique(round.(lon_all; digits = 4)))
ulat = sort(unique(round.(lat_all; digits = 4)))
ilon = [searchsortedfirst(ulon, round(x; digits = 4)) for x in lon]
ilat = [searchsortedfirst(ulat, round(x; digits = 4)) for x in lat]
function to_map(x)
    out = fill(NaN, length(ulon), length(ulat))
    for i in 1:nl
        ok[i] && (out[ilon[i], ilat[i]] = x[i])
    end
    return out
end
land = fill(NaN, length(ulon), length(ulat))
for i in 1:nl
    land[ilon[i], ilat[i]] = 1.0
end
ratio = [s1[i] > 0.5 ? log2(m1[i] / s1[i]) : NaN for i in 1:nl]
println("nodes with SoilGrids top metre below 0.5 kg C m⁻² (no ratio): ", count(i -> ok[i] && !(s1[i] > 0.5), 1:nl))

# zonal means over the land of each 1° latitude row
function zonal(x)
    out = fill(NaN, length(ulat))
    for j in eachindex(ulat)
        sel = [i for i in 1:nl if ok[i] && ilat[i] == j]
        isempty(sel) || (out[j] = sum(x[sel]) / length(sel))
    end
    return out
end

# global stock above each depth (Pg C)
depths = sort(unique(vcat(exp10.(range(log10(0.01), log10(15.0); length = 300)), [1.0, 2.0])))
okw = reshape(area .* ok, 1, :)
stock(x, D) = sum(replace(x, NaN => 0.0) .* above(D) .* okw) * 1e-12
stock_model = [stock(soc, D) for D in depths]
depths_sg = filter(<=(2.0), depths)
stock_sg = [stock(soilgrids, D) for D in depths_sg]

# figure
blues = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
diverging = ["#104281", "#2a78d6", "#86b6ef", "#f0efec", "#f2a3a2", "#e34948", "#9b2423"]
c_model, c_sg = "#2a78d6", "#eb6834"
ink, ink2, muted, grid_c = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1"
land_gray = "#c9c8c2"
fmt(x) = replace(@sprintf("%.0f", x), r"(\d)(?=(\d{3})+$)" => s"\1,")
g = totals["global"]

CM.set_theme!(fontsize = 15, textcolor = ink, Axis = (xgridcolor = grid_c, ygridcolor = grid_c, leftspinecolor = muted, bottomspinecolor = muted, topspinevisible = false, rightspinevisible = false, xtickcolor = muted, ytickcolor = muted, titlefont = :regular, titlealign = :left))
fig = CM.Figure(size = (1800, 1120), backgroundcolor = "#fcfcfb")
CM.Label(fig[0, 1:3], "ClimaLand soil organic carbon after the 1850 spin-up, against SoilGrids"; fontsize = 22, halign = :left, tellwidth = false)
CM.Label(fig[1, 1:3], "State: $note. Land north of 60°S. Top metre $(fmt(g[1])) vs $(fmt(g[2])) Pg C (×$(@sprintf("%.2f", g[1] / g[2]))); 0–2 m $(fmt(g[3])) vs $(fmt(g[4])) Pg C (×$(@sprintf("%.2f", g[3] / g[4]))); full 15 m column $(fmt(g[5])) Pg C, $(fmt(g[5] - g[3])) of it below 2 m, where SoilGrids has no data."; fontsize = 15, color = ink2, halign = :left, tellwidth = false, word_wrap = true)

map_axis(pos, title) = CM.Axis(pos; title, width = 560, height = 240, limits = ((-180, 180), (-60, 90)), xticks = -180:60:180, yticks = -60:30:90, xgridvisible = false, ygridvisible = false, xtickformat = v -> ["$(round(Int, x))°" for x in v], ytickformat = v -> ["$(round(Int, x))°" for x in v])
function map_panel!(pos, title, values, colors, colorrange)
    ax = map_axis(pos, title)
    colormap = CM.cgrad(colors)
    CM.heatmap!(ax, ulon, ulat, land; colormap = [land_gray, land_gray], colorrange = (0, 2))
    hm = CM.heatmap!(ax, ulon, ulat, values; colormap, colorrange, nan_color = :transparent, highclip = colors[end], lowclip = colors[1])
    return ax, hm
end
seq = blues
_, hm = map_panel!(fig[2, 1], "(a) ClimaLand, 0–1 m", to_map(m1), seq, (0, 40))
map_panel!(fig[2, 2], "(b) SoilGrids on the model grid, 0–1 m", to_map(s1), seq, (0, 40))
CM.Colorbar(fig[3, 1:2], hm; vertical = false, label = "SOC, 0–1 m (kg C m⁻²; gray: land left out)", width = CM.Relative(0.5), flipaxis = false)
_, hr = map_panel!(fig[2, 3], "(c) ClimaLand / SoilGrids, 0–1 m", to_map(ratio), diverging, (-2, 2))
CM.Colorbar(fig[3, 3], hr; vertical = false, label = "ratio (gray: SoilGrids below 0.5 kg C m⁻²)", ticks = (-2:2, ["1/4", "1/2", "1", "2", "4"]), width = CM.Relative(0.8), flipaxis = false)

ax = CM.Axis(fig[4, 1]; title = "(d) Zonal mean over land", xlabel = "latitude", ylabel = "SOC (kg C m⁻²)", xticks = -60:30:90, xtickformat = v -> ["$(round(Int, x))°" for x in v], limits = ((-60, 90), (0, 90)), width = 560, height = 300)
CM.lines!(ax, ulat, zonal(m1); color = c_model, linewidth = 2, label = "ClimaLand 0–1 m")
CM.lines!(ax, ulat, zonal(s1); color = c_sg, linewidth = 2, label = "SoilGrids 0–1 m")
CM.lines!(ax, ulat, zonal(m2); color = c_model, linewidth = 2, linestyle = :dash, label = "ClimaLand 0–2 m")
CM.lines!(ax, ulat, zonal(s2); color = c_sg, linewidth = 2, linestyle = :dash, label = "SoilGrids 0–2 m")
CM.axislegend(ax; position = :rt, framevisible = false, labelsize = 13, patchsize = (28, 10))
z2 = zonal(m2)
jmax = argmax(replace(z2, NaN => -Inf))
latname(x) = x < 0 ? "$(round(Int, -x))°S" : "$(round(Int, x))°N"
z2[jmax] > 90 && (ax.title = "(d) Zonal mean over land (ClimaLand 0–2 m clipped: $(round(Int, z2[jmax])) at $(latname(ulat[jmax])))")

deep = to_map(mall .- m2)
deep_max = sort(filter(isfinite, vec(deep)))[round(Int, 0.99 * count(isfinite, deep))]
_, hd = map_panel!(fig[4, 2], "(e) ClimaLand below 2 m (no SoilGrids data)", deep, seq, (0, 60))
println("99th percentile of SOC below 2 m: ", round(deep_max; digits = 1), " kg C m⁻²")
CM.Colorbar(fig[5, 2], hd; vertical = false, label = "SOC, 2–15 m (kg C m⁻²)", width = CM.Relative(0.8), flipaxis = false)

ax = CM.Axis(fig[4, 3]; title = "(f) Global stock above each depth", xlabel = "SOC above the depth (Pg C)", ylabel = "depth (m)", yscale = log10, yreversed = true, yticks = ([0.01, 0.1, 0.3, 1, 2, 5, 15], ["0.01", "0.1", "0.3", "1", "2", "5", "15"]), limits = ((0, nothing), (0.01, 15)), width = 560, height = 300)
CM.hlines!(ax, [1.0, 2.0]; color = muted, linewidth = 1, linestyle = :dot)
CM.lines!(ax, stock_model, depths; color = c_model, linewidth = 2, label = "ClimaLand")
CM.lines!(ax, stock_sg, depths_sg; color = c_sg, linewidth = 2, label = "SoilGrids (to 2 m)")
CM.text!(ax, stock_model[end], 15.0; text = fmt(stock_model[end]), align = (:right, :bottom), offset = (-4, 4), color = ink2, fontsize = 13)
CM.text!(ax, stock_sg[end], 2.0; text = fmt(stock_sg[end]), align = (:right, :top), offset = (-6, -4), color = ink2, fontsize = 13)
CM.axislegend(ax; position = :rt, framevisible = false, labelsize = 13, patchsize = (28, 10))
CM.rowgap!(fig.layout, 8)
CM.resize_to_layout!(fig)
CM.save(joinpath(out_dir, "soc_spinup_1850_vs_soilgrids.png"), fig; px_per_unit = 1.5)
println("wrote ", joinpath(out_dir, "soc_spinup_1850_vs_soilgrids.png"))

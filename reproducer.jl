import CairoMakie
import GeoMakie
import ClimaAnalysis

sim_var = OutputVar("global_diagnostics/output_0000/shf_1M_average.nc")
sim_var = ClimaAnalysis.average_time(sim_var)
sim_var = ClimaAnalysis.apply_oceanmask(sim_var)
data = sim_var.data
# data = randn(404, 202)

# struct DummyData{A}
#     data::A
# end
# dummy = DummyData(rand(404, 202))
# data = dummy.data

fig = CairoMakie.Figure(; size = (600, 400))
CairoMakie.Label(
        fig[1, 0],
        "shf",
        tellheight = false,
        fontsize = 30,
    )
layout = fig[1, 1] = CairoMakie.GridLayout()

lon = range(-180.0, 180.0, 404) |> collect
lat = range(-90.0, 90.0, 202) |> collect
title = "Sensible Heat Flux, average within 1 Month averaged over time\n(2.6784e6 to 5.995296e8s)"
ax = GeoMakie.GeoAxis(layout[1, 1]; title)
Makie.contourf!(ax, lon, lat, data)

CairoMakie.Label(fig[0, 1], "ANN", tellwidth = false, fontsize = 30)

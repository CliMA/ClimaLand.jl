using NCDatasets
using JLD2
using Plots

include("distributions.jl")
data = NCDataset(
    "/Users/katherinedeck/Downloads/6b0c4358-2bf3-4924-aa8f-793d468b92be(1)/ga2.nc",
)
lat = data["lat"][:];
lon = data["lon"][:];

sg_lat = (lat .< 34.4) .& (lat .> 34.16)
sg_lon = (lon .> -117.9) .& (lon .< -117.53)
Plots.heatmap(
    lon[sg_lon],
    lat[sg_lat],
    transpose(data["Band1"][sg_lon, sg_lat]),
    aspect_ratio = :equal,
)
mt_baldy = [-117.64641, 34.28885]
Plots.scatter!([mt_baldy[1]], [mt_baldy[2]], label = "Mt Baldy")
telegraph = [-117.59834, 34.26154]
Plots.scatter!([telegraph[1]], [telegraph[2]], label = "Telegraph Peak")
bp = [-117.76417, 34.35917]
Plots.scatter!([bp[1]], [bp[2]], label = "Baden-Powell")
sgr = [-117.85103, 34.21890]
Plots.scatter!([sgr[1]], [sgr[2]], label = "SG Res")
Plots.savefig("./topo_index_local.png")
ϕ = data["Band1"][sg_lon, sg_lat][:]
dist_type = InvGammaDistribution()
present = .~(typeof.(ϕ) .<: Missing) .& (ϕ .> 0)
present_count = sum(present)
x̄, var, params = fit(dist_type, ϕ[present])
ll_inv_gamma = sum(log.(pdf.(Ref(dist_type), ϕ, Ref(params))))
Plots.histogram(
    ϕ[present],
    label = "",
    xlabel = "Topographic index",
    ylabel = "Frequency",
    norm = true,
)
x = minimum(ϕ[present]):0.1:maximum(ϕ[present])
Plots.plot!(x, pdf.(Ref(dist_type), x, Ref(params)), label = "Inverse Gamma")
dist_type = GammaDistribution()
x̄, var, params = fit(dist_type, ϕ[present])
ll_gamma = sum(log.(pdf.(Ref(dist_type), ϕ, Ref(params))))
Plots.plot!(x, pdf.(Ref(dist_type), x, Ref(params)), label = "Gamma")
Plots.savefig("topo_local_distribution.png")

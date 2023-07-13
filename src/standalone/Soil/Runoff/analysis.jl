using JLD2
using Plots
include("distributions.jl")

files = ["./inv_gamma_1.0.jld2", "./lognormal_1.0.jld2", "./gamma_1.0.jld2","./frechet_1.0.jld2"]
p1 = load(files[1])["parameters"];
p2 = load(files[2])["parameters"];
p3 = load(files[3])["parameters"];
p4 = load(files[4])["parameters"];
LL1 = p1[:, :, 5]
LL2 = p2[:, :, 5]
LL3 = p3[:, :, 5]
LL4 = p4[:, :, 5]

Δ12 = LL1 .- LL2
Δ13 = LL1 .- LL3
Δ14 = LL1 .- LL4
plt = Plots.plot(xlabel = "ΔLL")
Plots.histogram!(plt, Δ12[Δ12 .!= 0.0], norm = true, label = "inv gamma/gamma")
Plots.histogram!(plt, Δ13[Δ13 .!= 0.0], norm = true, label = "inv gamma/lognormal")
Plots.histogram!(plt, Δ14[Δ14 .!= 0.0], norm = true, label = "inv gamma/frechet")
Plots.savefig(plt, "llΔ.png")

using NCDatasets
data = NCDataset(
    "/Users/katherinedeck/Downloads/6b0c4358-2bf3-4924-aa8f-793d468b92be(1)/ga2.nc",
)
lat = data["lat"][:];
lon = data["lon"][:];
(lat_min, lat_max) = extrema(lat)
(lon_min, lon_max) = extrema(lon)
(lat_min, lat_max) = extrema(lat)
(lon_min, lon_max) = extrema(lon)
resolution = 1.0
lat_count = Int(ceil((lat_max - lat_min) / resolution)) + 1
lon_count = Int(ceil((lon_max - lon_min) / resolution)) + 1
dist_type = InvGammaDistribution()
normalized = (LL1 ./ (eps(Float64) .+ p1[:, :, 6]))
c = findall(x -> x < -2.75, normalized)
for idc in c
    lat_id = idc[1]
    lon_id = idc[2]
    lat_indices =
        (lat .>= lat_min + resolution * (lat_id - 1)) .&
        (lat .< lat_min + resolution * lat_id)
    lon_indices =
        (lon .>= lon_min + resolution * (lon_id - 1)) .&
        (lon .< lon_min + resolution * lon_id)
    ϕ = data["Band1"][lon_indices, lat_indices][:]
    present = .~(typeof.(ϕ) .<: Missing) .& (ϕ .> 0)
    present_count = sum(present)
    Plots.histogram(ϕ[present], norm = true)
    α = p1[lat_id, lon_id, 8]
    β = p1[lat_id, lon_id, 9]
    Plots.scatter!(
        ϕ[present],
        pdf.(Ref(InvGammaDistribution()), ϕ[present], Ref([α, β])),
    )
    Plots.savefig()
end

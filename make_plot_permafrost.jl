using Statistics
using NCDatasets
using CairoMakie
root_path = "soil_longrun_gpu_varying_type/global_diagnostics/output_0000"
swc_data = NCDataset("$(root_path)/swc_1M_average.nc")
si_data = NCDataset("$(root_path)/si_1M_average.nc")
tsoil_data = NCDataset("$(root_path)/tsoil_1M_average.nc")
tair_data = NCDataset("$(root_path)/tair_1M_average.nc")
ssr_data = NCDataset("$(root_path)/ssr_1M_average.nc")
sr_data = NCDataset("$(root_path)/sr_1M_average.nc")
longlat = (16.5, 65.7)#(140.0, 58.0);#(-71.65, 52.23);#(-128.3, 69.1);#(-71.65, 52.23)#(143.2, 63.25)#
lat_id = findmin(abs.(swc_data["lat"][:] .- longlat[2]))[2]
lon_id = findmin(abs.(swc_data["lon"][:] .- longlat[1]))[2]


#snowc_data = NCDataset("$(root_path)/snowc_1M_average.nc")
#snd_data = NCDataset("$(root_path)/snd_1M_average.nc")
#iwc_data = NCDataset("$(root_path)/iwc_1M_average.nc")

times = swc_data["time"][:] ./(365*24*3600);
nelements = size(swc_data["swc"])[end]
fig = CairoMakie.Figure(size = (800, 1200))
ax1 = CairoMakie.Axis(fig[1, 1], xlabel = "Time", ylabel = "Tsoil - 273.15")
for i in 1:2:nelements
    lines!(
        ax1,
        tsoil_data["tsoil"][:,lon_id, lat_id, i] .- 273.15,
        label = "$i",
    )
end

ax2 = CairoMakie.Axis(fig[2, 1], xlabel = "Time", ylabel = "Ice")
for i in 1:2:nelements
    lines!(
        ax2,
        si_data["si"][:,lon_id, lat_id, i],
        label = "$i",
    )
end
fig[2, 2] = Legend(fig, ax2)

ax3 = CairoMakie.Axis(fig[3, 1], xlabel = "Time", ylabel = "Total water")
for i in 1:2:nelements
    lines!(
        ax3,
        si_data["si"][:,lon_id, lat_id, i] .+ swc_data["swc"][:,lon_id, lat_id, i],
        label = "$i",
    )
end
CairoMakie.save(joinpath(root_path, "soil_$(lon_id)_$(lat_id).png"), fig)
#=
fig = CairoMakie.Figure(size = (800, 1200))
ax1 = CairoMakie.Axis(fig[1, 1], xlabel = "Time", ylabel = "Snow depth")
lines!(
    ax1,
    snd_data["snd"][:,lon_id, lat_id],
)
ax2 = CairoMakie.Axis(fig[2, 1], xlabel = "Time", ylabel = "Snow cover")
lines!(
    ax2,
    snowc_data["snowc"][:,lon_id, lat_id],
)
ax3 = CairoMakie.Axis(fig[3, 1], xlabel = "Time", ylabel = "IWC")
lines!(
    ax3,
    iwc_data["iwc"][:,lon_id, lat_id],
)
CairoMakie.save(joinpath(root_path, "snow_$(lon_id)_$(lat_id).png"), fig)

fig = CairoMakie.Figure(size = (800, 1200))
ax1 = CairoMakie.Axis(fig[1, 1], xlabel = "Time", ylabel = "SSR")
lines!(
    ax1,
    ssr_data["ssr"][:,lon_id, lat_id],
)
ax2 = CairoMakie.Axis(fig[2, 1], xlabel = "Time", ylabel = "SR")
lines!(
    ax2,
    sr_data["sr"][:,lon_id, lat_id],
)
CairoMakie.save(joinpath(root_path, "runoff_$(lon_id)_$(lat_id).png"), fig)

last_two = collect((12*8+1):1:120)
T = tsoil_data["tsoil"][last_two,lon_id, lat_id,:];
ΔT_max = maximum(T, dims = 1) .- 273.15;
ΔT_max[ΔT_max .> 0] .= 1000
id = findmin(abs.(ΔT_max))[2][2]
depth = tsoil_data["z"][:][id]
T_min = findmin(T[:,id])[1]
σ = snowc_data["snowc"][last_two,lon_id, lat_id]
scdays = sum(σ .> 0.25)/length(σ)*365
@show depth
@show T_min
@show scdays
@show maximum(snd_data["snd"][last_two[σ .> 0], lon_id, lat_id])

=#

#active layer
times = sr_data["time"][:]
time_id = collect(length(times)-12-11:1:length(times))
T = tsoil_data["tsoil"][time_id,:, :,:];
T_air = tair_data["tair"][time_id,:, :];
annual_mean_tair = mean(T_air, dims = 1)[1,:,:] .- 273.15;
annual_mean_tsoil = mean(T, dims = 1)[1,:,:,:] .- 273.15;
annual_mean_swc = mean(swc_data["swc"][time_id,:,:,15], dims = 1)[1,:,:];
annual_mean_si = mean(si_data["si"][time_id,:,:,15], dims = 1)[1,:,:];

fig, ax, hm = heatmap(annual_mean_swc[:,110:end])
Colorbar(fig[:, end+1], hm)
scatter!(ax, [lon_id], [lat_id - 110], color=:black)
fig
save(joinpath(root_path, "annual_mean_swc.png"), fig)
fig, ax, hm = heatmap(annual_mean_si[:,110:end])
Colorbar(fig[:, end+1], hm)
scatter!(ax, [lon_id], [lat_id - 110], color=:black)
fig
save(joinpath(root_path, "annual_mean_si.png"), fig)

fig, ax, hm = heatmap(annual_mean_tair[:,110:end], colormap = :vik100, colorrange = (-30,30)) # maps of minimum ΔT at permafrost depth
Colorbar(fig[:, end+1], hm)
scatter!(ax, [lon_id], [lat_id - 110], color=:black)
fig
save(joinpath(root_path, "annual_mean_tair.png"), fig)

fig, ax, hm = heatmap((annual_mean_tsoil[:,:,15] .- annual_mean_tair)[:,110:end], colormap = :vik100, colorrange = (-5,5)) # maps of minimum ΔT at permafrost depth
Colorbar(fig[:, end+1], hm)
scatter!(ax, [lon_id], [lat_id - 110], color=:black)
fig
save(joinpath(root_path, "surface_deltaT.png"), fig)

fig, ax, hm = heatmap((annual_mean_tsoil[:,:,8] .- annual_mean_tair)[:,110:end], colorrange = (-5,5), colormap = :vik100) # maps of minimum ΔT at permafrost depth
Colorbar(fig[:, end+1], hm)
scatter!(ax, [lon_id], [lat_id - 110], color = :black)
fig
save(joinpath(root_path, "d2.5m_deltaT.png"), fig)
fig, ax, hm = heatmap((annual_mean_tsoil[:,:,4] .- annual_mean_tair)[:,110:end], colorrange = (-5,5), colormap = :vik100) # maps of minimum ΔT at permafrost depth
Colorbar(fig[:, end+1], hm)
scatter!(ax, [lon_id], [lat_id - 110], color = :black)
fig
save(joinpath(root_path, "d7m_deltaT.png"), fig)

fig, ax, hm = heatmap((annual_mean_tsoil[:,:,12] .- annual_mean_tair)[:,110:end], colorrange = (-5,5), colormap = :vik100) # maps of minimum ΔT at permafrost depth
Colorbar(fig[:, end+1], hm)
scatter!(ax, [lon_id], [lat_id - 110], color = :black)
fig
save(joinpath(root_path, "d12_deltaT.png"), fig)
fig, ax, hm = heatmap((annual_mean_tsoil[:,:,15] .- annual_mean_tair)[:,110:end], colorrange = (-5,5), colormap = :vik100) # maps of minimum ΔT at permafrost depth
Colorbar(fig[:, end+1], hm)
scatter!(ax, [lon_id], [lat_id - 110], color = :black)
fig
save(joinpath(root_path, "bottom_deltaT.png"), fig)


fig, ax, hm = heatmap((annual_mean_tsoil[:,:,15])[:,110:end], colormap = :vik100, colorrange = (-30,30)) # maps of minimum ΔT at permafrost depth
Colorbar(fig[:, end+1], hm)
scatter!(ax, [lon_id], [lat_id - 110], color = :black)
fig
save(joinpath(root_path, "surface.png"), fig)

fig, ax, hm = heatmap((annual_mean_tsoil[:,:,8])[:,110:end], colorrange = (-30,30), colormap = :vik100) # maps of minimum ΔT at permafrost depth
Colorbar(fig[:, end+1], hm)
scatter!(ax, [lon_id], [lat_id - 110], color = :black)
fig
save(joinpath(root_path, "d2.5m.png"), fig)
fig, ax, hm = heatmap((annual_mean_tsoil[:,:,4])[:,110:end], colorrange = (-30,30), colormap = :vik100) # maps of minimum ΔT at permafrost depth
Colorbar(fig[:, end+1], hm)
scatter!(ax, [lon_id], [lat_id - 110], color = :black)
fig
save(joinpath(root_path, "d7m.png"), fig)


ΔT_max = maximum(T, dims = 1) .- 273.15; # max at each layer and location - freezing. If this is < 0, that location is frozen for a full year at least one year
ΔT_max[ΔT_max .> 0] .= 1000; 
f(x) = findmin(x)[2]
ids = mapslices(f, abs.(ΔT_max), dims = 4)[1,:,:,:]; # find which layer each location maximum ΔT is at its minimum value - this is the layer closest to the freezing point and reflects the permafrost depth
cids = [CartesianIndex(i,j,ids[i,j,1]) for i in 1:360 for j in 1:180]
T_min = reshape(findmin(T[:,cids], dims = 1)[1], (180, 360))
T_min[T_min .> 273.15] .= NaN
fig, ax, hm = heatmap(transpose(T_min)[:,110:end], colorrange = (273.15-30,273.15)) # maps of minimum ΔT at permafrost depth
Colorbar(fig[:, end+1], hm)
fig
save(joinpath(root_path, "Tmin_north.png"), fig)
soil_depth = tsoil_data["z"][:]
ald = [soil_depth[ids[i,j,1]] for i in 1:360 for j in 1:180]
ald = reshape(ald, (180, 360))
fig2, ax2, hm2 = heatmap(transpose(ald)[:,110:end])
Colorbar(fig2[:, end+1], hm2)
fig2
save(joinpath(root_path, "depth_north.png"), fig2)
# Correlate permafrost depth and temperature with snow cover days
#=
σ = snowc_data["snowc"][:,:, :];
scdays = sum(σ .> 0.25, dims = 1)[1,:,:] ./size(σ)[1]*365
fig3, ax3, hm3 = heatmap(scdays[:,110:end])
Colorbar(fig3[:, end+1], hm3)
fig3
save(joinpath(root_path, "scdays_north.png"), fig3)
=#


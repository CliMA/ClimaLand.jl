using NCDatasets
using CairoMakie
using Statistics
root_path = "."
swc_data = NCDataset("$(root_path)/swc_1M_average.nc")
si_data = NCDataset("$(root_path)/si_1M_average.nc")
tsoil_data = NCDataset("$(root_path)/tsoil_1M_average.nc")
tair = NCDataset("$(root_path)/tair_1M_average.nc")

longlat =(-71.65, 50.23)# (-71.65, 52.23)#(-100.0,65.0)#(-71.65, 52.23)#(143.2, 63.25)#(-128.3, 69.1);
lat_id = findmin(abs.(swc_data["lat"][:] .- longlat[2]))[2]
lon_id = findmin(abs.(swc_data["lon"][:] .- longlat[1]))[2]
times = swc_data["time"][:] ./(365*24*3600);
nelements = size(swc_data["swc"])[end]
fig = CairoMakie.Figure(size = (800, 1200))
ax1 = CairoMakie.Axis(fig[1, 1], xlabel = "Time", ylabel = "T - 273.15")
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

fig = CairoMakie.Figure(size = (800, 800))
ax1 = CairoMakie.Axis(fig[1, 1], xlabel = "Time", ylabel = "T - mean(T_air)")
for i in 1:2:nelements
    lines!(
        ax1,
        tsoil_data["tsoil"][:,lon_id, lat_id, i] .- mean(tair["tair"][:,lon_id, lat_id]),
        label = "$i",
    )
end
CairoMakie.save(joinpath(root_path, "tsoil_rel_air_$(lon_id)_$(lat_id).png"), fig)



T = tsoil_data["tsoil"][240-12-11:end,:, :,:];
ΔT_max = maximum(T, dims = 1) .- 273.15; # max at each layer and location - freezing. If this is < 0, that location is frozen for a full year at least one year
ΔT_max[ΔT_max .> 0] .= 1000; 
f(x) = findmin(x)[2]
ids = mapslices(f, abs.(ΔT_max), dims = 4)[1,:,:,:]; # find which layer each location maximum ΔT is at its minimum value - this is the layer closest to the freezing point and reflects the permafrost depth
cids = [CartesianIndex(i,j,ids[i,j,1]) for i in 1:360 for j in 1:180]
T_min = reshape(findmin(T[:,cids], dims = 1)[1], (180, 360))
T_min[T_min .> 273.15] .= NaN
fig, ax, hm = heatmap(transpose(T_min)[:,110:end])
Colorbar(fig[1, 2], hm)
fig
save(joinpath(root_path, "tmin.png"), fig)

i = 9
mean_T_air = mean(tair["tair"][240-12*5-11:240,:,:], dims = 1)[1,:,:]
mean_T_sfc = mean(tsoil_data["tsoil"][240-12*5-11:240,:,:,i], dims = 1)[1,:,:]
mean_T_sfc[mean_T_air .>273.15] .=NaN
mean_T_air[mean_T_air .>273.15] .=NaN

lat = swc_data["lat"][:]
lon = swc_data["lon"][:]

fig = CairoMakie.Figure(size = (800,1600))
ax, hm = heatmap(fig[1,1], lon, lat[110:end], mean_T_air[:,110:end] .- 273.15,colorrange = (-30,0))
ax.limits = ((-180,180,lat[110], lat[end]))
ax.xticks = (-180:20:180, ["", "-160", "-140", "-120", "-100", "-80","-60", "-40", "-20", "0", "20","40","60","80","100","120","140","160",""])
ax.title= "Mean Air Temp- 273.15"
Colorbar(fig[1, 2], hm)
ax2, hm2 = heatmap(fig[2,1], lon, lat[110:end], mean_T_sfc[:,110:end] .- mean_T_air[:,110:end], colorrange=(-2,2), colormap=:vik) # maps of minimum ΔT at permafrost depth
Colorbar(fig[2, 2], hm2)
ax2.title= "mean T $(i) - mean T air"
ax2.limits = ((-180,180,lat[110], lat[end]))
ax2.xticks = (-180:20:180, ["", "-160", "-140", "-120", "-100", "-80","-60", "-40", "-20", "0", "20","40","60","80","100","120","140","160",""])

fig
save(joinpath(root_path, "T_$i.png"), fig)

m_w = 1000 .* sum(swc_data["swc"][:,:,:,:], dims = 4) .+ 917 .* sum(si_data["si"][:,:,:,:], dims = 4);

fig = CairoMakie.Figure(size = (800,1600))
ax, hm = heatmap(fig[1,1],sum(m_w[180-11:180,:,:,1], dims = 1)[1,:,:]./1000 ./15 .- sum(m_w[120-11:120,:,:,1], dims = 1)[1,:,:]./1000 ./15, colormap = :vik, colorrange = (-0.5,0.5))
ax.title= "∫θ(15 year) - ∫θ(10 year)"
Colorbar(fig[1, 2], hm)
ax2, hm2 = heatmap(fig[2,1], sum(m_w[240-11:240,:,:,1], dims = 1)[1,:,:]./1000 ./15 .- sum(m_w[180-11:180,:,:,1], dims = 1)[1,:,:]./1000 ./15, colormap=:vik, colorrange = (-0.5,0.5)) # maps of minimum ΔT at permafrost depth
Colorbar(fig[2, 2], hm2)
ax2.title= "∫θ(20 year) - ∫θ(15 year)"
fig

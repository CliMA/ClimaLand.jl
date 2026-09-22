using NCDatasets
using ArchGDAL
using SparseArrays
using JLD2

function make_basin_dict(;stem = "/Users/katherinedeck/Downloads/")
    pt_dataset = ArchGDAL.read(joinpath(stem,"hybas_pour_lev02_v1_shp/hybas_pour_lev02_v1.shp"))
    basin_dict = Dict()
    for f in ArchGDAL.getlayer(pt_dataset, 0)
        id =  Int(ArchGDAL.getfield(f, 0))
        geom = ArchGDAL.getgeom(f)
        x = round(ArchGDAL.getx(geom, 0))
        y = round(ArchGDAL.gety(geom, 0))
        coords = (x,y)
        item= (outlets = [ArchGDAL.createpoint(coords...),],)
        if haskey(basin_dict, id)
            if ArchGDAL.createpoint(coords...) ∈ basin_dict[id].outlets
                nothing
            else
                push!(basin_dict[id].outlets, ArchGDAL.createpoint(coords...))
            end
        else
            basin_dict[id] = item
        end
    end

    prefixes = ["af","ar", "as", "au", "eu", "si", "na", "sa", "gr"]
    for prefix in prefixes
        @show prefix
        dataset = ArchGDAL.read(joinpath(stem, "hybas_$(prefix)_lev02_v1c/hybas_$(prefix)_lev02_v1c.shp"))
        layer = ArchGDAL.getlayer(dataset, 0)
        for feature in layer
            id = Int(ArchGDAL.getfield(feature, 0))
            geom = ArchGDAL.getgeom(feature)
            nt = (; geom = geom, endo = ArchGDAL.getfield(feature,9), area = ArchGDAL.getfield(feature, 6), centroid = ArchGDAL.centroid(geom))
            basin_dict[id] = merge(basin_dict[id], nt)
        end
    end
    # Antartica does not have a basin...
    antartica_id = 1
    antartica_geom = ArchGDAL.createpolygon([(-181.0, -60.0), (181.0, -60.0), (180.0, -91.0), (-181.0, -91.0), (-181.0, -60.0)])
    basin_dict[antartica_id] = (;geom = antartica_geom, endo = 1, area = 14.2e6, outlets = [ArchGDAL.createpoint((0.0,-90.0))], centroid = ArchGDAL.centroid(antartica_geom))
    return basin_dict
end
function check_containment(x,y; geom)
    pt = ArchGDAL.createpoint(x,y); return ArchGDAL.contains(geom, pt)
end


function make_grid_to_basin_map(x,y, basin_dict, pixel_area)
    npixels = length(x)
    tmp = zeros(length(x))
    rowvals = []
    nzvals = []
    colptrs = []
    nbasins = length(keys(basin_dict))
    R_earth = 6371000.0/1000.0# in km
    for (idx,key) in enumerate(keys(basin_dict))
        @show idx
        geom = basin_dict[key].geom
        ch(x,y) = check_containment(x,y;geom = geom)
        tmp .= ch.(x,y)
        grid_area = @. tmp* pixel_area*cos(y*π/180) *(π/180)^2 * R_earth^2 # instead do sum over CC field, which will compute area correctly? But N.B that for the global box domain this is in degrees ^2. How to handle missed pixels then?
        basin_area = basin_dict[key].area
        N = Int(sum(tmp))
        push!(rowvals, findall(tmp[:] .==1))
        push!(colptrs, zeros(Int, N) .+ idx)
        push!(nzvals, grid_area[tmp .==1]./basin_area)
        tmp .= 0
    end
    mapped = vcat(rowvals...);
    columns =  vcat(colptrs...,)
    values = vcat(nzvals...,)
    regrid_matrix = sparse(mapped, columns,values, npixels, nbasins);
    # foo = ones(npixels)
    #∑ f A = ∑f_b A_b 
    # sum(@. 1 * pixel_area* cos(y*π/180) *(π/180)^2 * R_earth^2) ≈sum((reshape(foo, (1,22420)) * regrid_matrix)[:] .* basin_areas) # true
    return regrid_matrix
end


P = NCDataset("precip_1M_average.nc")
ET  =NCDataset("et_1M_average.nc")
TR  =NCDataset("tr_1M_average.nc")
TWSA(P,ET,TR) = P+ET+TR
lon = P["lon"][:]
lat = P["lat"][:]

tmp = P["precip"][1,:,:]
mask = .~isnan.(tmp)
x = reshape(repeat(lon, outer = 180), (360,180))[mask]
y = reshape(repeat(lat, inner = 360), (360,180))[mask]
basin_dict = make_basin_dict()
regrid_matrix_1deg = make_grid_to_basin_map(x,y,basin_dict,1.0);

twsa = NCDataset("twsa_0.5x0.5.nc")
grace_mask = .~ismissing.(twsa["twsa"][:,:,4]);
grace_x  = reshape(repeat(twsa["lon"], outer = 360), (720,360))[grace_mask];
grace_y = reshape(repeat(twsa["lat"], inner = 720), (720,360))[grace_mask];
regrid_matrix_halfdeg = make_grid_to_basin_map(grace_x,grace_y,basin_dict,0.5*0.5);
@save "twsa_comparison_nomissing.jld2" regrid_matrix_1deg regrid_matrix_halfdeg basin_dict
# save to JLD2 file; takes forever to make regrid matrix

cl = (reshape(ones(22420), (1,22420)) * regrid_matrix_1deg)[:] .* basin_areas
gr = (reshape(ones(72652), (1,72652)) * regrid_matrix_halfdeg)[:] .* basin_areas
area_mismatch = findall(abs.(gr .- cl) ./ cl .> 0.1)
keep = [k for k in 1:63 if ~(k∈findall(abs.(gr .- cl) ./ cl .> 0.1))]
clima_basin = zeros(length(keep), 156)
grace_basin = zeros(length(keep), 156)
for i in 1:156
    k = i+22
    clima = -1 .* 30*86400*(TR["tr"][k,:,:] .* 1000 .+ P["precip"][k,:,:] .+ ET["et"][k, :,:]); # seconds per month not general)
    clima_basin[:,i] .= (reshape(clima[mask], (1, 22420)) * regrid_matrix_1deg)[:][keep]
   # @show Second(P["time"][:][i+22]) .+ Month(2) .+ DateTime(2000)
   # @show twsa["time"][:][i]
    if length(unique(twsa["twsa"][:,:,i])) > 1
        grace = twsa["twsa"][:,:,i];
        grace_basin[:,i] .= (reshape(grace[grace_mask], (1, 72652)) * regrid_matrix_halfdeg)[:][keep]
    else grace_basin[:,i] .= NaN
    end
end
basin_ids = [k for k in keys(basin_dict)]
nanmean(x) = mean(x[.~isnan.(x)])
nanmedian(x) = median(x[.~isnan.(x)])

amp_grace = mapslices(nanmedian, abs.(grace_basin), dims = 2).- mapslices(nanmean, grace_basin, dims = 2)
amp_clima = mapslices(nanmedian, abs.(clima_basin), dims = 2) .- mapslices(nanmean, clima_basin, dims = 2)
RMSE = sqrt.(mapslices(nanmean, (clima_basin .- grace_basin).^2, dims = 2))

fig = Figure(); ax = Axis(fig[1,1])
for (idx,id) in enumerate(basin_ids[keep])
    plot!(ax, basin_dict[id].geom, color = min(Int(round(RMSE[idx])),100), colormap = :viridis, colorrange = (0,100))
end
Colorbar(fig[1,2], limits = (0,100), colormap = :viridis)

fig = Figure(size= (1500,800)); ax = Axis(fig[1,1], title = "Clima");ax2 = Axis(fig[1,2], title = "GRACE")
for (idx,id) in enumerate(basin_ids[keep])
    plot!(ax, basin_dict[id].geom, color = min(Int(round(amp_clima[idx])),150), colormap = :viridis, colorrange = (0,150))
    plot!(ax2, basin_dict[id].geom, color = min(Int(round(amp_grace[idx])),150), colormap = :viridis, colorrange = (0,150))
end


Colorbar(fig[1,3], limits = (0,150), colormap = :viridis)


fig = Figure(size= (1500,800)); ax = Axis(fig[1,1], title = "Clima");ax2 = Axis(fig[1,2], title = "GRACE")
for (idx,id) in enumerate(basin_ids)
    plot!(ax, basin_dict[id].geom, color = Int(round(cl[idx]/1e5)), colormap = :viridis, colorrange = (0,100))
    plot!(ax2, basin_dict[id].geom, color = Int(round(gr[idx]/1e5)), colormap = :viridis, colorrange = (0,100))
end


Colorbar(fig[1,3], limits = (0,100), colormap = :viridis)


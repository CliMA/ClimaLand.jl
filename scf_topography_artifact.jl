using NCDatasets
using Statistics
inpath = "/Users/katherinedeck/.julia/artifacts/fe19d8dbe7a18ff39588e1f718014b0479d9c0f7/ETOPO_2022_v1_60s_N90W180_surface.nc"
data = NCDataset(inpath)

lat = data["lat"][:];
dlat = lat[2] - lat[1]
lon = data["lon"][:];
dlon = lon[2] - lon[1]
z = data["z"][:,:];

function slope(lat, dlat,dlon, dzx, dzy)
    dx = dlat*π/180 * R_earth # dx refers to distance spanned by dlat
    dy = dlon*π/180 * cos(lat*π/180) * R_earth # dy refers to distnace spanned by dlon
    μs = 1/2*((dzx/dx)^2  + (dzy/dy)^2) # squared slope, unitless sinc z is in meters
    return μs
end

R_earth =6371000 # in meters
μs = zeros(21600,10800);
for j in 1:21600
    @show j
    for i in 1:10800
        # Get dz in `y` direction (longitude)
        if j == 1 
            dzy = z[j, i] - z[end, i] # periodic in longitude
        else
            dzy = z[j, i] - z[j-1, i]
        end
        # Get dz in "x" direction (latitude)
        if i == 1
            dzx = z[j,i+1] - z[j, i] # approx using  value at i = 2
        else
            dzx = z[j,i] - z[j, i-1]
        end

        μs[j, i] = slope(lat[i], dlat,dlon,  dzx, dzy)
    end
end

# Now compute the statistics we need for bins of size res x res
# μ = sqrt(mean(μs))
# σz = std(z)
# ξ = σz/μ/L
# mean(z)
# Bin to 1/2 degree
res = 0.5 # 360/res must be an integer; 180/res must be an integer
outlat = (-90.0+res/2):res:(90.0-res/2)
outlong = (-180.0+res/2):res:(180.0-res/2)
lon_mask = BitArray(undef, 21600)
lat_mask = BitArray(undef, 10800)
nlon = Int(360/res)
nlat = Int(180/res)
ξ = zeros(nlon, nlat)
μ = zeros(nlon, nlat)
mean_z = zeros(nlon, nlat)
for (i,δ) in enumerate(outlat)
    for (j,λ) in enumerate(outlong)
        lon_mask .= (lon .< (λ + res/2)) .&& (lon .>= (λ -res/2))
        lat_mask .= (lat .< (δ + res/2)) .&& (lat .>= (δ -res/2))
        masked_z = z[lon_mask, lat_mask]
        σ_z = std(masked_z)
        dA = R_earth^2*(π/180)^2*res^2(1+ cos(δ)^2)
        L = sqrt(dA)
        μ[j,i] = sqrt(mean(μs[lon_mask, lat_mask]))
        ξ[j,i] = σ_z/L/max(μ[j,i], eps(Float32))
        mean_z[j,i] = mean(masked_z)
    end
end
close(data)

outpath = "topographic_data.nc"
ds = NCDataset(outpath,"c")
defDim(ds,"lon",nlon)
defDim(ds,"lat",nlat)
la = defVar(ds, "lat", Float32, ("lat",))
lo = defVar(ds, "lon", Float32, ("lon",))
la.attrib["units"] = "degrees_north"
la.attrib["standard_name"] = "latitude"
lo.attrib["standard_name"] = "longitude"
lo.attrib["units"] = "degrees_east"
la[:] = outlat
lo[:] = outlong
# Define a global attribute
ds.attrib["title"] = "Topographic Statistics for subgrid parameterizations"

# Define the variables temperature
v1 = defVar(ds,"mean_height",Float32,("lon","lat"))
v1[:,:] = mean_z
v2 = defVar(ds,"norm_stdz",Float32,("lon","lat"))
v2[:,:] = ξ
v3 = defVar(ds,"mean_slope",Float32,("lon","lat"))
μmax = quantile(μ[:], 0.9999) # remove extrema
v3[:,:] = min.(μ, μmax)

v1.attrib["units"] = "meters"
v1.attrib["comments"] = "Average height of surface in grid cell"
v2.attrib["units"] = ""
v2.attrib["comments"] = "Standard deviation of height normalized by slope and cell size; measure of complexity of topography inside the cell"
v3.attrib["units"] = ""
v3.attrib["comments"] = "Average slope of surface in grid cell"
close(ds)

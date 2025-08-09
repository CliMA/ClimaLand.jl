# temporary file for retrieving site-specific fluxnet params 


"""
    get_site_parameters(site_ID::String)

Returns a named tuple containing all site-specific parameters for the given site.
Currently supports: "US-Var", "US-Ha1", "US-NR1", "US-MOz"

Only includes parameters that are actually used in the forward_model_site function.
"""
function get_site_parameters(site_ID::String)
    if site_ID == "US-Var"
        return (
            # Site information
            data_link = "https://caltech.box.com/shared/static/54huhy74kxnn23i6w1s54vddo2j4hl71.csv",
            time_offset = 8,
            lat = FT(38.4133),
            long = FT(-120.9508),
            atmos_h = FT(2),
            
            # Domain parameters
            n_stem = Int64(0),
            n_leaf = Int64(1), 
            h_leaf = FT(0.5),
            h_stem = FT(0),
            nelements = 14,
            zmin = FT(-0.5),
            zmax = FT(0),
            dz_tuple = nothing,
            N_days = 364,
            
            # Soil parameters
            soil_ν = FT(0.45),
            soil_K_sat = FT(0.45 / 3600 / 100),
            soil_S_s = FT(1e-3),
            soil_vg_n = FT(2.0),
            soil_vg_α = FT(2.0),
            θ_r = FT(0.067),
            ν_ss_quartz = FT(0.3),
            ν_ss_om = FT(0.02),
            ν_ss_gravel = FT(0.0),
            z_0m_soil = FT(0.01),
            z_0b_soil = FT(0.01),
            soil_ϵ = FT(0.98),
            soil_α_PAR = FT(0.2),
            soil_α_NIR = FT(0.2),
            
            # TwoStreamModel parameters
            Ω = FT(1.0),
            ld = FT(0.5),
            G_Function = ConstantGFunction(FT(0.5)),
            α_PAR_leaf = FT(0.11),
            λ_γ_PAR = FT(5e-7),
            τ_PAR_leaf = FT(0.05),
            α_NIR_leaf = FT(0.35),
            τ_NIR_leaf = FT(0.34),
            ϵ_canopy = FT(0.97),
            
            # Conductance Model
            g1 = FT(166),
            Drel = FT(1.6),
            g0 = FT(1e-4),
            
            # Photosynthesis model
            is_c3 = FT(1.0),
            Vcmax25 = FT(2 * 4.225e-5),
            
            # Energy Balance model
            ac_canopy = FT(745),
            
            # Plant Hydraulics parameters
            pc = FT(-3e5),
            sc = FT(1e-3),
            SAI = FT(0),
            f_root_to_shoot = FT(1.0),
            K_sat_plant = 2e-8,
            ψ63 = FT(-2.7 / 0.0098),
            Weibull_param = FT(4),
            a = FT(0.05 * 0.0098),
            conductivity_model = PlantHydraulics.Weibull{FT}(2e-8, FT(-2.7 / 0.0098), FT(4)),
            retention_model = PlantHydraulics.LinearRetentionCurve{FT}(FT(0.05 * 0.0098)),
            plant_ν = FT(8.93e-3),
            plant_S_s = FT(1e-2 * 0.0098),
            rooting_depth = FT(0.3),
            
            # Aerodynamic parameters
            h_canopy = max(FT(0.5), FT(0)),
            z0_m = FT(0.13) * max(FT(0.5), FT(0)),
            z0_b = FT(0.1) * FT(0.13) * max(FT(0.5), FT(0)),
        )
    elseif site_ID == "US-Ha1"
        return (
            # Site information
            data_link = "https://caltech.box.com/shared/static/xixaod6511cutz51ag81k1mtvy05hbol.csv",
            time_offset = 5,
            lat = FT(42.5378),
            long = FT(-72.1715),
            atmos_h = FT(30),
            
            # Domain parameters
            n_stem = Int64(1),
            n_leaf = Int64(1),
            h_leaf = FT(12),
            h_stem = FT(14),
            nelements = 20,
            zmin = FT(-10),
            zmax = FT(0),
            dz_tuple = (FT(1.5), FT(0.025)),
            N_days = 364,
            
            # Soil parameters
            soil_ν = FT(0.5),
            soil_K_sat = FT(4e-7),
            soil_S_s = FT(1e-3),
            soil_vg_n = FT(2.05),
            soil_vg_α = FT(0.04),
            θ_r = FT(0.067),
            ν_ss_quartz = FT(0.1), 
            ν_ss_om = FT(0.1),
            ν_ss_gravel = FT(0.0),
            z_0m_soil = FT(0.01),
            z_0b_soil = FT(0.001),
            soil_ϵ = FT(0.98),
            soil_α_PAR = FT(0.2),
            soil_α_NIR = FT(0.2),
            
            # TwoStreamModel parameters
            Ω = FT(0.69),
            ld = FT(0.5),
            G_Function = ConstantGFunction(FT(0.5)),
            α_PAR_leaf = FT(0.1),
            λ_γ_PAR = FT(5e-7),
            τ_PAR_leaf = FT(0.05),
            α_NIR_leaf = FT(0.45),
            τ_NIR_leaf = FT(0.25),
            ϵ_canopy = FT(0.97),
            
            # Conductance Model
            g1 = FT(141),
            Drel = FT(1.6),
            g0 = FT(1e-4),
            
            # Photosynthesis model
            is_c3 = FT(1.0),
            Vcmax25 = FT(9e-5),
            
            # Energy Balance model
            ac_canopy = FT(2.5e3),
            
            # Plant Hydraulics parameters
            SAI = FT(1.0),
            f_root_to_shoot = FT(3.5),
            K_sat_plant = 5e-9,
            ψ63 = FT(-4 / 0.0098),
            Weibull_param = FT(4),
            a = FT(0.05 * 0.0098),
            conductivity_model = PlantHydraulics.Weibull{FT}(5e-9, FT(-4 / 0.0098), FT(4)),
            retention_model = PlantHydraulics.LinearRetentionCurve{FT}(FT(0.05 * 0.0098)),
            plant_ν = FT(2.46e-4),
            plant_S_s = FT(1e-2 * 0.0098),
            rooting_depth = FT(0.5),
        )
    elseif site_ID == "US-NR1"
        return (
            # Site information
            data_link = "https://caltech.box.com/shared/static/r6gvldgabk3mvtx53gevnlaq1ztsk41i.csv",
            time_offset = 7,
            lat = FT(40.0329),
            long = FT(-105.5464),
            atmos_h = FT(21.5),
            
            # Domain parameters
            n_stem = Int64(1),
            n_leaf = Int64(1),
            h_leaf = FT(6.5),
            h_stem = FT(7.5),
            nelements = 20,
            zmin = FT(-10),
            zmax = FT(0),
            dz_tuple = (FT(1.25), FT(0.05)),
            N_days = 364,
            
            # Soil parameters
            soil_ν = FT(0.45),
            soil_K_sat = FT(4e-7),
            soil_S_s = FT(1e-3),
            soil_vg_n = FT(2.05),
            soil_vg_α = FT(0.04),
            θ_r = FT(0.0),
            ν_ss_quartz = FT(0.1),
            ν_ss_om = FT(0.1),
            ν_ss_gravel = FT(0.0),
            z_0m_soil = FT(0.1),
            z_0b_soil = FT(0.1),
            soil_ϵ = FT(0.98),
            soil_α_PAR = FT(0.2),
            soil_α_NIR = FT(0.2),
            
            # TwoStreamModel parameters
            Ω = FT(0.71),
            ld = FT(0.5),
            G_Function = ConstantGFunction(FT(0.5)),
            α_PAR_leaf = FT(0.1),
            λ_γ_PAR = FT(5e-7),
            τ_PAR_leaf = FT(0.05),
            α_NIR_leaf = FT(0.35),
            τ_NIR_leaf = FT(0.25),
            ϵ_canopy = FT(0.97),
            
            # Conductance Model
            g1 = FT(141),
            Drel = FT(1.6),
            g0 = FT(1e-4),
            
            # Photosynthesis model
            is_c3 = FT(1.0),
            Vcmax25 = FT(9e-5),
            
            # Energy Balance model
            ac_canopy = FT(3e3),
            
            # Plant Hydraulics parameters
            SAI = FT(1.0),
            f_root_to_shoot = FT(3.5),
            K_sat_plant = 5e-9,
            ψ63 = FT(-4 / 0.0098),
            Weibull_param = FT(4),
            a = FT(0.05 * 0.0098),
            conductivity_model = PlantHydraulics.Weibull{FT}(5e-9, FT(-4 / 0.0098), FT(4)),
            retention_model = PlantHydraulics.LinearRetentionCurve{FT}(FT(0.05 * 0.0098)),
            plant_ν = FT(8.06e-4),
            plant_S_s = FT(1e-2 * 0.0098),
            rooting_depth = FT(1.0),
        )
    elseif site_ID == "US-MOz"
        return (
            # Site information
            data_link = "https://caltech.box.com/shared/static/7r0ci9pacsnwyo0o9c25mhhcjhsu6d72.csv",
            time_offset = 7,
            lat = FT(38.7441),
            long = FT(-92.2000),
            atmos_h = FT(32),
            
            # Domain parameters
            n_stem = Int64(1),
            n_leaf = Int64(1),
            h_leaf = FT(9.5),
            h_stem = FT(9),
            nelements = 20,
            zmin = FT(-10),
            zmax = FT(0),
            dz_tuple = (FT(1.5), FT(0.1)),
            N_days = 364,
            
            # Soil parameters
            soil_ν = FT(0.55),
            soil_K_sat = FT(4e-7),
            soil_S_s = FT(1e-2),
            soil_vg_n = FT(2.0),
            soil_vg_α = FT(0.05),
            θ_r = FT(0.04),
            ν_ss_quartz = FT(0.1),
            ν_ss_om = FT(0.1),
            ν_ss_gravel = FT(0.0),
            z_0m_soil = FT(0.01),
            z_0b_soil = FT(0.01),
            soil_ϵ = FT(0.98),
            soil_α_PAR = FT(0.2),
            soil_α_NIR = FT(0.2),
            
            # TwoStreamModel parameters
            Ω = FT(0.69),
            χl = FT(0.1),
            G_Function = CLMGFunction(FT(0.1)),
            α_PAR_leaf = FT(0.1),
            λ_γ_PAR = FT(5e-7),
            τ_PAR_leaf = FT(0.05),
            α_NIR_leaf = FT(0.45),
            τ_NIR_leaf = FT(0.25),
            ϵ_canopy = FT(0.97),
            
            # Conductance Model
            g1 = FT(141),
            Drel = FT(1.6),
            g0 = FT(1e-4),
            
            # Photosynthesis model
            is_c3 = FT(1.0),
            Vcmax25 = FT(6e-5),
            
            # Energy Balance model
            ac_canopy = FT(5e2),
            
            # Plant Hydraulics parameters
            pc = FT(-2.0e6),
            sc = FT(5e-6),
            SAI = FT(1.0),
            f_root_to_shoot = FT(3.5),
            K_sat_plant = 7e-8,
            ψ63 = FT(-4 / 0.0098),
            Weibull_param = FT(4),
            a = FT(0.1 * 0.0098),
            conductivity_model = PlantHydraulics.Weibull{FT}(7e-8, FT(-4 / 0.0098), FT(4)),
            retention_model = PlantHydraulics.LinearRetentionCurve{FT}(FT(0.1 * 0.0098)),
            plant_ν = FT(1.44e-4),
            plant_S_s = FT(1e-2 * 0.0098),
            rooting_depth = FT(0.5),
        )
    else
        error("Site ID $site_ID not supported. Available sites: US-Var, US-Ha1, US-NR1, US-MOz")
    end
end


using NCDatasets

"""
    nearest_point_value(path, varname, qlat, qlon, qdepth)

Return the value of `varname` at the grid point nearest to the given
latitude (`qlat`), longitude (`qlon`), and depth (`qdepth`) from the NetCDF
file at `path`.

Assumes:
- Rectilinear grid (1D lat, 1D lon, 1D depth)
- No time dimension
- Variable includes latitude, longitude, and depth dimensions in any order

Only the coordinate vectors and the single selected value are read from disk.
"""
function nearest_point_value(path::AbstractString, varname::AbstractString,
                             qlat::Real, qlon::Real, qdepth::Real)

    to360(x) = (x % 360 + 360) % 360
    to180(x) = x > 180 ? x - 360 : x
    nearest_index(v, q) = argmin(abs.(v .- q))

    ds = NCDataset(path, "r")
    try
        # Identify dimension positions by common aliases
        v = ds[varname]
        dims = dimnames(v)  # e.g., ("depth","lat","lon")
        latdim = findfirst(in(("lat","latitude","y","nav_lat")), dims)
        londim = findfirst(in(("lon","longitude","x","nav_lon")), dims)
        depthdim = findfirst(in(("depth","lev","level","z","olevel","Depth")), dims)

        latdim === nothing && error("Could not identify latitude dimension in $varname.")
        londim === nothing && error("Could not identify longitude dimension in $varname.")
        depthdim === nothing && error("Could not identify depth/level dimension in $varname.")

        # Coordinate variable names usually match dim names in CF; use them directly
        lat = ds[dims[latdim]][:]
        lon = ds[dims[londim]][:]
        depth = ds[dims[depthdim]][:]

        # Adjust query longitude to file’s convention
        lonq = maximum(lon) > 180 ? to360(qlon) : to180(qlon)

        i_lat = nearest_index(lat, qlat)
        i_lon = nearest_index(lon, lonq)
        i_depth = nearest_index(depth, qdepth)

        idx = Any[Colon() for _ in 1:length(dims)]
        idx[latdim] = i_lat
        idx[londim] = i_lon
        idx[depthdim] = i_depth

        return v[idx...]
    finally
        close(ds)
    end
end
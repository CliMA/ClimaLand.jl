module Rivers
using ClimaLand
using ArchGDAL
using ClimaCore
using StaticArrays
using DocStringExtensions

import ClimaLand:
    AbstractExpModel,
    name,
    auxiliary_vars,
    auxiliary_types,
    auxiliary_domain_names,
    make_update_boundary_fluxes

using ClimaLand.Domains
import ClimaLand.Domains:
    coordinates
export RiverModel

abstract type AbstractRiverModel{FT} <: AbstractExpModel{FT} end

"""
    InstantaneousRouting{FT, M, D, RD} <: AbstractRiverModel{FT}

$(DocStringExtensions.FIELDS)
"""
struct InstantaneousRouting{FT, GD, BD, F} <: AbstractRiverModel{FT}
    "The grid for the land model"
    grid::GD
    "The river basin dictionary"
    basin_dict::BD
    "Matrix summing the gridded flux to basin: x_grid * M = x_basin"
    regrid_matrix::F
end

function InstantaneousRouting{FT}(cc_domain, basin_dict) where {FT}
    mask = parent(cc_domain.space.surface.grid.mask.is_active)[:]
    coords = ClimaLand.Domains.coordinates(cc_domain).surface
    y = parent(coords.:1)[:][mask]
    x = parent(coords.:2)[:][mask]
    grid = (;x,y, mask)
    regrid_matrix = make_grid_to_basin_map(grid, basin_dict)
    args = (grid, basin_dict, regrid_matrix)
    return InstantaneousRouting{FT, typeof.(args)...}(args...)
end

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
function make_grid_to_basin_map(grid, basin_dict)
    (;x,y) = grid
    npixels = length(x)
    tmp = zeros(length(x))
    rowvals = []
    nzvals = []
    colptrs = []
    nbasins = length(keys(basin_dict))
    pixel_area = 1.0 # need to generalize - how?
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
    # Missing some points along coastlines, inland seas (Caspian, Black Sea...) since land points don't align with basins boundaries exactly.
    possible =  collect(1:1:npixels)
    missed = [s for s in possible if ~(s∈ mapped)]
    missed_rows = []
    missed_cols = []
    missed_vals = []
    basin_ids = [k for k in keys(basin_dict)]
    basin_areas = [basin_dict[k].area for k in keys(basin_dict)]
    for pixel in missed
        basin_distance_list = []
        pt = ArchGDAL.createpoint(x[pixel], y[pixel])
        for key in basin_ids
            centroid = basin_dict[key].centroid
            dist = ArchGDAL.distance(pt, centroid) # in degrees
            push!(basin_distance_list, dist)
        end
        id = findmin(basin_distance_list)
        push!(missed_cols, id[2])
        basin_area = basin_dict[basin_ids[id[2]]].area
        push!(missed_vals, pixel_area  *(π/180)^2 * R_earth^2 *  cos(y[pixel]*π/180)/basin_area)
    end
    
    all_indices = vcat(mapped..., missed...)
    columns =  vcat(colptrs..., missed_cols...)
    values = vcat(nzvals...,missed_vals...)
    regrid_matrix = sparse(all_indices, columns,values, npixels, nbasins);
    # foo = ones(npixels)
    # sum(@. cos(y*π/180) *(π/180)^2 * R_earth^2) ≈sum((reshape(foo, (1,22420)) * regrid_matrix)[:] .* basin_areas) # true
    return regrid_matrix
end

ClimaLand.name(model::InstantaneousRouting) = :river

ClimaLand.auxiliary_vars(model::InstantaneousRouting) = (:R_basin,:R_tot_grid)
ClimaLand.auxiliary_types(model::InstantaneousRouting{FT}) where {FT} = (FT,FT)
ClimaLand.auxiliary_domain_names(model::InstantaneousRouting) = (:basin, :grid)

function Domain.coordinates(basin_dict, FT)
    return @SVector zeros(FT, length(basin_dict))
end

function Domain.coordinates(grid, FT)
    N = sum(grid.mask)
    return @SVector zeros(FT, N)
end

function Domains.coordinates(model::InstantaneousRouting{FT}) where {FT}
    npixels = sum(model.grid.mask)
    pixel_coords = @SVector zeros(FT, npixels)
    nbasins = length(model.basin_dict)
    basin_coords = @SVector zeros(FT, nbasins)
    return (;grid = pixel_coords, basin = basin_coords)# this breaks - why?
end

function make_update_boundary_fluxes(model::InstantaneousRouting)
    function update_boundary_fluxes(p,Y,t)
        R_tot_grid = p.rivers.R_tot_grid 
        R_tot_grid .= parent(p.soil.R_s)[:][mask] .+ parent(p.soil.R_ss)[:][mask]
        p.rivers.R_basin .= reshape(R_tot_grid, (1,length(R_tot_grid)) * regrid_matrix # Takes 0.00006 seconds for 1 degree land
    end
    return update_boundary_fluxes
end

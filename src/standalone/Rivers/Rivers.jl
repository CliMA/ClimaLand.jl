module Rivers
using ClimaLand
using StaticArrays
using ArchGDAL
using DocStringExtensions

import ClimaLand:
    AbstractExpModel,
    name,
    auxiliary_vars,
    auxiliary_types,
    auxiliary_domain_names,
    make_update_boundary_fluxes,

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
    "The ClimaCore grid (domain) for the land model"
    grid::GD
    "The river basin dictionary"
    basin_dict::BD
    "Matrix summing the gridded flux to basin: x_grid * M = x_basin"
    regrid_matrix::F
    function InstantaneousRouting{FT}(grid_domain, basin_domain) where {FT}
        regrid_matrix = make_grid_to_basin_map(sum, grid_domain, basin_domain)
        args = (grid_domain, basin_domain, regrid_matrix)
        return new{FT, typeof.(args)...}(args...)
    end
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
            nt = (; geom = ArchGDAL.getgeom(feature), endo = ArchGDAL.getfield(feature,9), area = ArchGDAL.getfield(feature, 6))
            basin_dict[id] = merge(basin_dict[id], nt)
        end
    end
    return basin_dict
end

function make_grid_to_basin_map(grid, basin_dict)
    coords = ClimaLand.Domains.coordinates(grid).surface
#    y = parent(coords.:1)[:]
#    x = parent(coords.:2)[:]
#    mask = parent(axes(coords).grid.mask.is_active)[:]
#    y = y[mask]
#    x = x[mask]
    #    tmp = zero(x)
    y = coords.:1
    x = coords.:2
    tmp = ClimaCore.Fields.zeros(x)
    rowvals = []
    nzvals = []
    colptrs = []
    nbasins = length(keys(basin_dict))
    for (idx,key) in enumerate(keys(basin_dict))
        @show idx
        geom = basin_dict[key].geom
        ch(x,y) = check_containment(x,y;geom = geom)
        tmp .= ch.(x,y)
        
        grid_area = sum(tmp)
        basin_area = basin_dict[key].area
        N = Int(sum(parent(tmp)))
        push!(rowvals, findall(parent(tmp)[:] .==1)) # full length, 60k entries. Not masked to land
        push!(colptrs, zeros(Int, N) .+ idx)
        push!(nzvals, zeros(N).+grid_area/basin_area)
        tmp .= 0
    end
    mapped = vcat(rowvals...);
    # Missing Antartica, some points along coastlines, inland seas (Caspian, Black Sea...)
    # Currently just map to a fake "basin"
    possible =  collect(1:1:length(x))
    missed = [s for s in posible if ~(s∈ mapped)]
    all_indices = (mappped..., missed...)
    columns =  vcat(colptrs..., ones(Int, length(missed)) .+ nbasins)
    values = vcat(nzvals..., ones(length(missed))) # this is not correct...
    regrid_matrix = sparse(all_indices, columnsd,values, length(x), nbasins+1);
    
    return regrid_matrix
end

ClimaLand.name(model::InstantaneousRouting) = :river

ClimaLand.auxiliary_vars(model::InstantaneousRouting) = (:R_outlet,:R_tot_grid)
ClimaLand.auxiliary_types(model::InstantaneousRouting{FT}) where {FT} = (FT,FT)
ClimaLand.auxiliary_domain_names(model::InstantaneousRouting) = (:basin, :grid)

function Domain.coordinates(basin_domain::SVector)
    return river_coordinates(basin_domain)
end

function Domain.coordinates(model::InstantaneousRouting)
    return (:grid = coordinates(model.grid_domain), :basin = coordinates(model.basin_domain))
end

function make_update_boundary_fluxes(model::InstantaneousRouting)
    function update_boundary_fluxes(p,Y,t)
        @. p.rivers.R_tot_grid = p.soil.R_s + p.soil.R_ss
        tmp = parent(p.rivers.R_tot_grid)[mask] # allocates, need to define mask...
        p.rivers.R_outlet .= reshape(x, (1,length(x))) * regrid_matrix # Takes 0.00006 seconds for 1 degree land
    end
    return update_boundary_fluxes
end

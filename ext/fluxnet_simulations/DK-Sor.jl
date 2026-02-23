#=
This file contains domain specifications, parameters, and other site-level information for running
ClimaLand at the DK-Sor (Soroe, Denmark) fluxtower site.
DK-Sor is a temperate deciduous broadleaf (beech) forest, FLUXNET2015 Tier 1 site.
Citation: Pilegaard, K., et al. (2011), Increasing net CO2 uptake by a Danish beech forest during
the period from 1996 to 2009, Agricultural and Forest Meteorology, 151(7), 934-946.
=#

"""
    get_domain_info(FT, ::Val{:DK_Sor}; kwargs...)

Gets and returns primary domain information for the DK-Sor (Soroe, Denmark) site,
which is a temperate deciduous broadleaf (beech) forest.
"""
function FluxnetSimulations.get_domain_info(
    FT,
    ::Val{:DK_Sor};
    dz_bottom = FT(2.0),
    dz_top = FT(0.05),
    nelements = 20,
    zmin = FT(-10.0),
    zmax = FT(0),
)
    dz_tuple = (dz_bottom, dz_top)
    return (; dz_tuple, nelements, zmin, zmax)
end

"""
    get_location(FT, ::Val{:DK_Sor}; kwargs...)

Returns geographical information for DK-Sor (Soroe, Denmark) site.
The `time_offset` is the difference from UTC in hours
and excludes daylight savings time, following Fluxnet convention.
For this site, the local time is UTC+1 (CET).
"""
function FluxnetSimulations.get_location(
    FT,
    ::Val{:DK_Sor};
    time_offset = 1,
    lat = FT(55.48587),
    long = FT(11.64464),
)
    return (; time_offset, lat, long)
end

"""
    get_fluxtower_height(FT, ::Val{:DK_Sor}; kwargs...)

Returns atmosphere height for DK-Sor site.
Reference height = 57 m, canopy height = 25 m (from NetCDF metadata).
"""
function FluxnetSimulations.get_fluxtower_height(
    FT,
    ::Val{:DK_Sor};
    atmos_h = FT(57),
)
    return (; atmos_h,)
end


###################################
#        MODULE FUNCTIONS         #
###################################

"""
    get_sensor_depths(FT, ::Val{site_ID})

Return the depths (m) of the soil temperature (`tsoil`) and soil water content
(`swc`) sensors of the FLUXNET site, ordered as the `TS_F_MDS_i` and
`SWC_F_MDS_i` columns, as a NamedTuple of tuples; `nothing` for a sensor set
whose depths are not documented. Sites without a method fall back to unknown
depths.
"""
FluxnetSimulations.get_sensor_depths(FT, ::Val) =
    (; tsoil = nothing, swc = nothing)

###################################
#            UTILITIES            #
###################################

"""
    replace_hyphen(old_site_ID::String)

Replaces all instances of hyphens in a given site ID string with underscores
and returns a Symbol of the reformatted site ID to be used as a Val{} type.

For example, an input string "US-MOz" would be output as "US_MOz".
"""
function FluxnetSimulations.replace_hyphen(old_site_ID::String)
    new_site_ID = replace(old_site_ID, "-" => "_")

    return Symbol(new_site_ID)
end

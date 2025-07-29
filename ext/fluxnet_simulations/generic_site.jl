###################################
#        MODULE FUNCTIONS         #
###################################

# TODO

###################################
#            UTILITIES            #
###################################

"""
    replace_hyphen(old_site_ID::String)

Replaces all instances of hyphens in a given site ID string with underscores
and returns a Symbol of the reformatted site ID to be used as a Val{} type.

For example, an input string "US-MOz" would be output as "US_MOz".
"""
function replace_hyphen(old_site_ID::String)
    new_site_ID = replace(old_site_ID, "-" => "_")

    return Symbol(new_site_ID)
end

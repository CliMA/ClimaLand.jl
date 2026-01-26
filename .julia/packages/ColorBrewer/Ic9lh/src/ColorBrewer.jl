module ColorBrewer

using Colors
using JSON

# The colorbrewer.json file contains the schemes,
# and was obtained from colorbrewer2.org.
# The contents of that file are covered under the Apache License 2.0, included
# in this distribution under data/.
# Here, dirname is used to accommodate non-standard package directory
# configurations (as well as standard ones, obviously).
colorSchemes = JSON.parsefile(
    dirname(@__FILE__) * "/../data/colorbrewer.json");

# Creates a color palette based on the colorbrewer2 color schemes.
#
# Args:
#   schemeName: The name of the scheme. The schemes are as follows:
#       - sequential: OrRd, PuBu, BuPu, Oranges, BuGn, YlOrBr, YlGn, Reds,
#           RdPu, Greens, YlGnBu, Purples, GnBu, Greys, YlOrRd, PuRd, Blues,
#           PuBuGn
#       - diverging: Spectral, RdYlGn, RdBu, PiYG, PRGn, RdYlBu, BrBG, RdGy,
#           PuOr
#       - qualitative: Set1, Set2, Set3, Accent, Dark2, Paired, Pastel1, Pastel2
#   n: The number of colors in the palette. Must be within the proper range.
#       The range for each type is as follows:
#       - sequential: between 3 and 9, inclusive
#       - diverging: between 3 and 11, inclusive
#       - qualitative: between 3 and 12 for Paired, Set3
#                      between 3 and 9 for Pastel1, Set1
#                      between 3 and 8 for Set2, Accent, Dark2, Pastel2
#
# Throws:
#   ArgumentError: If the scheme name is invalid.
#   ArgumentError: If the number of colors in the palette is invalid.
#
# Returns: An array of RGBs (from the Colors package).
function palette(schemeName::AbstractString, n::Integer)
    if !haskey(colorSchemes, schemeName)
        throw(ArgumentError("Scheme $schemeName is undefined."));
    end #Close if statement validating scheme name.

    if !haskey(colorSchemes[schemeName], string(n))
        throw(ArgumentError(
            "Scheme $schemeName is undefined for $n colors."));
    end #Close if statement validating n for the scheme name.

    return [parse(Colorant, s) for s in colorSchemes[schemeName][string(n)]];
end

export palette

end # module

module Readme
using Jieko: @pub, DEF
using DocStringExtensions: SIGNATURES, TYPEDSIGNATURES

const MyFancyAlias = Real

"""
$DEF

my lovely interface
"""
@pub jieko(x::MyFancyAlias) = x

"""
$SIGNATURES

my lovely method
"""
doc_string_ext(x::MyFancyAlias) = x


"""
$TYPEDSIGNATURES

my lovely method
"""
doc_string_ext_type(x::MyFancyAlias) = x

end # Readme

module Multi
using Jieko: @pub, DEF

const MyAliasName = Int

"""
$DEF
"""
@pub jieko(x::Real)::Int = error("not implemented")

"""
$DEF
"""
@pub jieko(x::MyAliasName)::Complex = error("not implemented")

end

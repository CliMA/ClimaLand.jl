"""
$DEFLIST
"""
module Basic

using Jieko: @pub, DEF, DEFLIST, @prelude_module

"""
$DEF

This is a constant.
"""
@pub const X = 1

"""
$DEF

This is a type.
"""
@pub struct Foo <: Real
    x::Int
end

"""
$DEF

This is a macro.
"""
@pub macro goo(x::Int, y::String, zs...)
end

"""
$DEF

This is a macro.
"""
@pub macro moo(x::Int, y::String="aaa")
end

"""
$DEF

This is a function.
"""
@pub foo(x::Float64)::Int = 2

@prelude_module

end # module

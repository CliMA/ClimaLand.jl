# This file is a part of InverseFunctions.jl, licensed under the MIT License (MIT).
"""
    InverseFunctions

Lightweight package that defines an interface to invert functions.
"""
module InverseFunctions

include("functions.jl")
include("inverse.jl")
include("setinverse.jl")

"""
    InverseFunctions.test_inverse(f, x; compare=isapprox, kwargs...)

Test if [`inverse(f)`](@ref) is implemented correctly.

The function tests (as a `Test.@testset`) if

* `compare(inverse(f)(f(x)), x) == true` and
* `compare(inverse(inverse(f))(x), f(x)) == true`.

`kwargs...` are forwarded to `compare`.

!!! Note
    On Julia >= 1.9, you have to load the `Test` standard library to be able to use
    this function.
"""
function test_inverse end

@static if !isdefined(Base, :get_extension)
    include("../ext/InverseFunctionsDatesExt.jl")
    include("../ext/InverseFunctionsTestExt.jl") 
end

# Better error message if users forget to load Test
if isdefined(Base, :get_extension) && isdefined(Base.Experimental, :register_error_hint)
    function __init__()
        Base.Experimental.register_error_hint(MethodError) do io, exc, _, _
            if exc.f === test_inverse &&
                (Base.get_extension(InverseFunctions, :InverseFunctionsTest) === nothing)
                print(io, "\nDid you forget to load Test?")
            end
        end
    end
end

end # module

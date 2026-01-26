# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Rmath_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Rmath")
JLLWrappers.@generate_main_file("Rmath", UUID("f50d1b31-88e8-58de-be2c-1cc44531875f"))
end  # module Rmath_jll

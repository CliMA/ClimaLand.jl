# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Imath_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Imath")
JLLWrappers.@generate_main_file("Imath", UUID("905a6f67-0a94-5f89-b386-d35d92009cd1"))
end  # module Imath_jll

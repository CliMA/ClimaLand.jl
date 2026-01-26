# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule LZO_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("LZO")
JLLWrappers.@generate_main_file("LZO", UUID("dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"))
end  # module LZO_jll

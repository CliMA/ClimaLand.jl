# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule OpenBLAS32_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("OpenBLAS32")
JLLWrappers.@generate_main_file("OpenBLAS32", UUID("656ef2d0-ae68-5445-9ca0-591084a874a2"))
end  # module OpenBLAS32_jll

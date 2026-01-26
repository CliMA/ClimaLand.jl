# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule NVTX_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("NVTX")
JLLWrappers.@generate_main_file("NVTX", UUID("e98f9f5b-d649-5603-91fd-7774390e6439"))
end  # module NVTX_jll

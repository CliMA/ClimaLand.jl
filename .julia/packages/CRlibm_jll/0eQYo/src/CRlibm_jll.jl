# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule CRlibm_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("CRlibm")
JLLWrappers.@generate_main_file("CRlibm", UUID("4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"))
end  # module CRlibm_jll

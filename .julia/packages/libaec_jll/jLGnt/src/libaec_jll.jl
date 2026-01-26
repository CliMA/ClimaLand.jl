# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule libaec_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("libaec")
JLLWrappers.@generate_main_file("libaec", UUID("477f73a3-ac25-53e9-8cc3-50b2fa2566f0"))
end  # module libaec_jll

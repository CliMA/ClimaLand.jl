# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule libfdk_aac_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("libfdk_aac")
JLLWrappers.@generate_main_file("libfdk_aac", UUID("f638f0a6-7fb0-5443-88ba-1cc74229b280"))
end  # module libfdk_aac_jll

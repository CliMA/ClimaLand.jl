# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule oneTBB_jll
using Base
using Base: UUID
using LazyArtifacts
import JLLWrappers

JLLWrappers.@generate_main_file_header("oneTBB")
JLLWrappers.@generate_main_file("oneTBB", UUID("1317d2d5-d96f-522e-a858-c73665f53c3e"))
end  # module oneTBB_jll

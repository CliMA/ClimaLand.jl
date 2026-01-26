# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule XZ_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("XZ")
JLLWrappers.@generate_main_file("XZ", UUID("ffd25f8a-64ca-5728-b0f7-c24cf3aae800"))
end  # module XZ_jll

# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule LERC_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("LERC")
JLLWrappers.@generate_main_file("LERC", UUID("88015f11-f218-50d7-93a8-a6af411a945d"))
end  # module LERC_jll

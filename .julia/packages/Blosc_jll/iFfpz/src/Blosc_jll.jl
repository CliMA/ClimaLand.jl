# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Blosc_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Blosc")
JLLWrappers.@generate_main_file("Blosc", UUID("0b7ba130-8d10-5ba8-a3d6-c5182647fed9"))
end  # module Blosc_jll

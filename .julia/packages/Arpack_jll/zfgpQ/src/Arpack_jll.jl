# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Arpack_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Arpack")
JLLWrappers.@generate_main_file("Arpack", UUID("68821587-b530-5797-8361-c406ea357684"))
end  # module Arpack_jll

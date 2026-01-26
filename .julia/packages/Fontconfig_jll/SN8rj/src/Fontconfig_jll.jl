# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Fontconfig_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Fontconfig")
JLLWrappers.@generate_main_file("Fontconfig", UUID("a3f928ae-7b40-5064-980b-68af3947d34b"))
end  # module Fontconfig_jll

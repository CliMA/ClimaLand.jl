# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Xorg_libXau_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Xorg_libXau")
JLLWrappers.@generate_main_file("Xorg_libXau", UUID("0c0b7dd1-d40b-584c-a123-a41640f87eec"))
end  # module Xorg_libXau_jll

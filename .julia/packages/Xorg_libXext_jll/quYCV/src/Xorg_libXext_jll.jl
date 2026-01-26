# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Xorg_libXext_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Xorg_libXext")
JLLWrappers.@generate_main_file("Xorg_libXext", UUID("1082639a-0dae-5f34-9b06-72781eeb8cb3"))
end  # module Xorg_libXext_jll

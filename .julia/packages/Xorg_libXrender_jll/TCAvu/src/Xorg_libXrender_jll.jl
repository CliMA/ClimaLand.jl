# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Xorg_libXrender_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Xorg_libXrender")
JLLWrappers.@generate_main_file("Xorg_libXrender", UUID("ea2f1a96-1ddc-540d-b46f-429655e07cfa"))
end  # module Xorg_libXrender_jll

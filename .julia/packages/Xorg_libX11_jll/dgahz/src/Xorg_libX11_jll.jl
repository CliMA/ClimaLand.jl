# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Xorg_libX11_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Xorg_libX11")
JLLWrappers.@generate_main_file("Xorg_libX11", UUID("4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"))
end  # module Xorg_libX11_jll

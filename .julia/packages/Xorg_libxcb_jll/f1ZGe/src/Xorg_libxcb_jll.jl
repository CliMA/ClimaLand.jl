# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Xorg_libxcb_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Xorg_libxcb")
JLLWrappers.@generate_main_file("Xorg_libxcb", UUID("c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"))
end  # module Xorg_libxcb_jll

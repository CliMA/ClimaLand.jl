# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Xorg_libXdmcp_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Xorg_libXdmcp")
JLLWrappers.@generate_main_file("Xorg_libXdmcp", UUID("a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"))
end  # module Xorg_libXdmcp_jll

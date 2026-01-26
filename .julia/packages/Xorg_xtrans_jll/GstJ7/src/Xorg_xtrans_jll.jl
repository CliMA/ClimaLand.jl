# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Xorg_xtrans_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Xorg_xtrans")
JLLWrappers.@generate_main_file("Xorg_xtrans", UUID("c5fb5394-a638-5e4d-96e5-b29de1b5cf10"))
end  # module Xorg_xtrans_jll

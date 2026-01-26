# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Xorg_libpciaccess_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Xorg_libpciaccess")
JLLWrappers.@generate_main_file("Xorg_libpciaccess", UUID("a65dc6b1-eb27-53a1-bb3e-dea574b5389e"))
end  # module Xorg_libpciaccess_jll

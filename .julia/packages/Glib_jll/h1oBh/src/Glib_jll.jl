# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Glib_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Glib")
JLLWrappers.@generate_main_file("Glib", UUID("7746bdde-850d-59dc-9ae8-88ece973131d"))
end  # module Glib_jll

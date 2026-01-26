# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Ogg_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Ogg")
JLLWrappers.@generate_main_file("Ogg", UUID("e7412a2a-1a6e-54c0-be00-318e2571c051"))
end  # module Ogg_jll

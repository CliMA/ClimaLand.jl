# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Graphite2_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Graphite2")
JLLWrappers.@generate_main_file("Graphite2", UUID("3b182d85-2403-5c21-9c21-1e1f0cc25472"))
end  # module Graphite2_jll

# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Expat_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Expat")
JLLWrappers.@generate_main_file("Expat", UUID("2e619515-83b5-522b-bb60-26c02a35a201"))
end  # module Expat_jll

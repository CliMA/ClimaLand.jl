# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Pango_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Pango")
JLLWrappers.@generate_main_file("Pango", UUID("36c8627f-9965-5494-a995-c6b170f724f3"))
end  # module Pango_jll

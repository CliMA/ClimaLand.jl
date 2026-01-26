# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule GettextRuntime_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("GettextRuntime")
JLLWrappers.@generate_main_file("GettextRuntime", UUID("b0724c58-0f36-5564-988d-3bb0596ebc4a"))
end  # module GettextRuntime_jll

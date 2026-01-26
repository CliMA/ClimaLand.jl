# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Hwloc_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Hwloc")
JLLWrappers.@generate_main_file("Hwloc", UUID("e33a78d0-f292-5ffc-b300-72abe9b543c8"))
end  # module Hwloc_jll

# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Libtiff_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Libtiff")
JLLWrappers.@generate_main_file("Libtiff", UUID("89763e89-9b03-5906-acba-b20f662cd828"))
end  # module Libtiff_jll

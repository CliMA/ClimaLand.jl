# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule x264_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("x264")
JLLWrappers.@generate_main_file("x264", UUID("1270edf5-f2f9-52d2-97e9-ab00b5d0237a"))
end  # module x264_jll

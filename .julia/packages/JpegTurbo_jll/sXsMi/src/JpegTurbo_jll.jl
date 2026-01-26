# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule JpegTurbo_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("JpegTurbo")
JLLWrappers.@generate_main_file("JpegTurbo", UUID("aacddb02-875f-59d6-b918-886e6ef4fbf8"))
end  # module JpegTurbo_jll

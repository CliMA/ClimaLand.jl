# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule isoband_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("isoband")
JLLWrappers.@generate_main_file("isoband", UUID("9a68df92-36a6-505f-a73e-abb412b6bfb4"))
end  # module isoband_jll

# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule OpenEXR_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("OpenEXR")
JLLWrappers.@generate_main_file("OpenEXR", UUID("18a262bb-aa17-5467-a713-aee519bc75cb"))
end  # module OpenEXR_jll

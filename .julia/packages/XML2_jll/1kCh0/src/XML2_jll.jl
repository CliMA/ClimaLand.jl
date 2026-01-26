# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule XML2_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("XML2")
JLLWrappers.@generate_main_file("XML2", UUID("02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"))
end  # module XML2_jll

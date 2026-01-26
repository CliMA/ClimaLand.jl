# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Cairo_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Cairo")
JLLWrappers.@generate_main_file("Cairo", UUID("83423d85-b0ee-5818-9007-b63ccbeb887a"))
end  # module Cairo_jll

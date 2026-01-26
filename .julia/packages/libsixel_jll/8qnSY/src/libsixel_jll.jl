# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule libsixel_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("libsixel")
JLLWrappers.@generate_main_file("libsixel", UUID("075b6546-f08a-558a-be8f-8157d0f608a5"))
end  # module libsixel_jll

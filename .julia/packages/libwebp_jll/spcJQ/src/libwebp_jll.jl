# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule libwebp_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("libwebp")
JLLWrappers.@generate_main_file("libwebp", UUID("c5f90fcd-3b7e-5836-afba-fc50a0988cb2"))
end  # module libwebp_jll

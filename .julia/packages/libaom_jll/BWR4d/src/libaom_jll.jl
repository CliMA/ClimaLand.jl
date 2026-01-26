# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule libaom_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("libaom")
JLLWrappers.@generate_main_file("libaom", UUID("a4ae2306-e953-59d6-aa16-d00cac43593b"))
end  # module libaom_jll

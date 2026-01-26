# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule LLVMOpenMP_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("LLVMOpenMP")
JLLWrappers.@generate_main_file("LLVMOpenMP", UUID("1d63c593-3942-5779-bab2-d838dc0a180e"))
end  # module LLVMOpenMP_jll

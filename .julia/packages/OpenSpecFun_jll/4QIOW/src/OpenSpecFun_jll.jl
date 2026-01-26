# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule OpenSpecFun_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("OpenSpecFun")
JLLWrappers.@generate_main_file("OpenSpecFun", UUID("efe28fd5-8261-553b-a9e1-b2916fc3738e"))
end  # module OpenSpecFun_jll

# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule EarCut_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("EarCut")
JLLWrappers.@generate_main_file("EarCut", UUID("5ae413db-bbd1-5e63-b57d-d24a61df00f5"))
end  # module EarCut_jll

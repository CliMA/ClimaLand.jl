# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Libmount_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Libmount")
JLLWrappers.@generate_main_file("Libmount", UUID("4b2f31a3-9ecc-558c-b454-b3730dcb73e9"))
end  # module Libmount_jll

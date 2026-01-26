# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Lz4_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Lz4")
JLLWrappers.@generate_main_file("Lz4", UUID("5ced341a-0733-55b8-9ab6-a4889d929147"))
end  # module Lz4_jll

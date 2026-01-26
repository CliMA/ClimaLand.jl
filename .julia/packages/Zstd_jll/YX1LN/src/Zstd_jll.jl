# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Zstd_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Zstd")
JLLWrappers.@generate_main_file("Zstd", UUID("3161d3a3-bdf6-5164-811a-617609db77b4"))
end  # module Zstd_jll

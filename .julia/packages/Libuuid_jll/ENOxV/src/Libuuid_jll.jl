# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Libuuid_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Libuuid")
JLLWrappers.@generate_main_file("Libuuid", UUID("38a345b3-de98-5d2b-a5d3-14cd9215e700"))
end  # module Libuuid_jll

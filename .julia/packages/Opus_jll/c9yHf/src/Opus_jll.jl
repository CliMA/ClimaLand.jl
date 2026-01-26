# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Opus_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Opus")
JLLWrappers.@generate_main_file("Opus", UUID("91d4177d-7536-5919-b921-800302f37372"))
end  # module Opus_jll

# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Giflib_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Giflib")
JLLWrappers.@generate_main_file("Giflib", UUID("59f7168a-df46-5410-90c8-f2779963d0ec"))
end  # module Giflib_jll

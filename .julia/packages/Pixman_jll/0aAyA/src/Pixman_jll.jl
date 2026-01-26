# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Pixman_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Pixman")
JLLWrappers.@generate_main_file("Pixman", UUID("30392449-352a-5448-841d-b1acce4e97dc"))
end  # module Pixman_jll

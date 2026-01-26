# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Libglvnd_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Libglvnd")
JLLWrappers.@generate_main_file("Libglvnd", UUID("7e76a0d4-f3c7-5321-8279-8d96eeed0f29"))
end  # module Libglvnd_jll

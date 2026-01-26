# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule LAME_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("LAME")
JLLWrappers.@generate_main_file("LAME", UUID("c1c5ebd0-6772-5130-a774-d5fcae4a789d"))
end  # module LAME_jll

# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule FreeType2_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("FreeType2")
JLLWrappers.@generate_main_file("FreeType2", UUID("d7e528f0-a631-5988-bf34-fe36492bcfd7"))
end  # module FreeType2_jll

# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule FFMPEG_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("FFMPEG")
JLLWrappers.@generate_main_file("FFMPEG", UUID("b22a6f82-2f65-5046-a5b2-351ab43fb4e5"))
end  # module FFMPEG_jll

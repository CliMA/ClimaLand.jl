# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule libvorbis_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("libvorbis")
JLLWrappers.@generate_main_file("libvorbis", UUID("f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"))
end  # module libvorbis_jll

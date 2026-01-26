# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Libiconv_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Libiconv")
JLLWrappers.@generate_main_file("Libiconv", UUID("94ce4f54-9a6c-5748-9c1c-f9c7231a4531"))
end  # module Libiconv_jll

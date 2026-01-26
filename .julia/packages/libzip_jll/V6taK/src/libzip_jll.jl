# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule libzip_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("libzip")
JLLWrappers.@generate_main_file("libzip", UUID("337d8026-41b4-5cde-a456-74a10e5b31d1"))
end  # module libzip_jll

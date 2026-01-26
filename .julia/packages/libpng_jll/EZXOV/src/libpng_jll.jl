# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule libpng_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("libpng")
JLLWrappers.@generate_main_file("libpng", UUID("b53b4c65-9356-5827-b1ea-8c7a1a84506f"))
end  # module libpng_jll

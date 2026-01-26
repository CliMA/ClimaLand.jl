# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule libass_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("libass")
JLLWrappers.@generate_main_file("libass", UUID("0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"))
end  # module libass_jll

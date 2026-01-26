# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule HarfBuzz_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("HarfBuzz")
JLLWrappers.@generate_main_file("HarfBuzz", UUID("2e76f6c2-a576-52d4-95c1-20adfe4de566"))
end  # module HarfBuzz_jll

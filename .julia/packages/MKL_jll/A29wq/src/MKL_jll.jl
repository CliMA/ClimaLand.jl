# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule MKL_jll
using Base
using Base: UUID
using LazyArtifacts
import JLLWrappers

JLLWrappers.@generate_main_file_header("MKL")
JLLWrappers.@generate_main_file("MKL", UUID("856f044c-d86e-5d09-b602-aeab76dc8ba7"))
end  # module MKL_jll

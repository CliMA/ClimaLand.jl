# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule IntelOpenMP_jll
using Base
using Base: UUID
using LazyArtifacts
import JLLWrappers

JLLWrappers.@generate_main_file_header("IntelOpenMP")
JLLWrappers.@generate_main_file("IntelOpenMP", UUID("1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"))
end  # module IntelOpenMP_jll

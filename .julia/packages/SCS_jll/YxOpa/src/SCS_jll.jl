# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule SCS_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("SCS")
JLLWrappers.@generate_main_file("SCS", UUID("f4f2fc5b-1d94-523c-97ea-2ab488bedf4b"))
end  # module SCS_jll

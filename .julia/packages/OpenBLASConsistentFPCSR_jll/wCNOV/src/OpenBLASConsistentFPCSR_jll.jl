# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule OpenBLASConsistentFPCSR_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("OpenBLASConsistentFPCSR")
JLLWrappers.@generate_main_file("OpenBLASConsistentFPCSR", UUID("6cdc7f73-28fd-5e50-80fb-958a8875b1af"))
end  # module OpenBLASConsistentFPCSR_jll

# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule OpenSSL_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("OpenSSL")
JLLWrappers.@generate_main_file("OpenSSL", UUID("458c3c95-2e84-50aa-8efc-19380b2a3a95"))
end  # module OpenSSL_jll

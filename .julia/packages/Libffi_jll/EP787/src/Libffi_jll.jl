# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Libffi_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Libffi")
JLLWrappers.@generate_main_file("Libffi", UUID("e9f186c6-92d2-5b65-8a66-fee21dc1b490"))
end  # module Libffi_jll

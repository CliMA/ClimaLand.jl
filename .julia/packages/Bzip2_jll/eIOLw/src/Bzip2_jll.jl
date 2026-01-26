# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule Bzip2_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("Bzip2")
JLLWrappers.@generate_main_file("Bzip2", UUID("6e34b625-4abd-537c-b88f-471c36dfa7a0"))
end  # module Bzip2_jll

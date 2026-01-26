# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule JuliaNVTXCallbacks_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("JuliaNVTXCallbacks")
JLLWrappers.@generate_main_file("JuliaNVTXCallbacks", UUID("9c1d0b0a-7046-5b2e-a33f-ea22f176ac7e"))
end  # module JuliaNVTXCallbacks_jll

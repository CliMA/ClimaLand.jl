# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule FriBidi_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("FriBidi")
JLLWrappers.@generate_main_file("FriBidi", UUID("559328eb-81f9-559d-9380-de523a88c83c"))
end  # module FriBidi_jll

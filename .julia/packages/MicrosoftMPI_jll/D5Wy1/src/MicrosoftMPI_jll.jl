# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule MicrosoftMPI_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("MicrosoftMPI")
JLLWrappers.@generate_main_file("MicrosoftMPI", UUID("9237b28f-5490-5468-be7b-bb81f5f5e6cf"))
end  # module MicrosoftMPI_jll

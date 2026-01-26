# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule FFTW_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("FFTW")
JLLWrappers.@generate_main_file("FFTW", UUID("f5851436-0d7a-5f13-b9de-f02708fd171a"))
end  # module FFTW_jll

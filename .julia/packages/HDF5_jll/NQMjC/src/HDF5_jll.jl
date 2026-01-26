# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule HDF5_jll
using Base
using Base: UUID
using LazyArtifacts
using MPIPreferences
Base.include(@__MODULE__, joinpath("..", ".pkg", "platform_augmentation.jl"))
import JLLWrappers

JLLWrappers.@generate_main_file_header("HDF5")
JLLWrappers.@generate_main_file("HDF5", UUID("0234f1f7-429e-5d53-9886-15a909be8d59"))
end  # module HDF5_jll

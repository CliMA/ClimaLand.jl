# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule NetCDF_jll
using Base
using Base: UUID
using LazyArtifacts
using MPIPreferences
Base.include(@__MODULE__, joinpath("..", ".pkg", "platform_augmentation.jl"))
import JLLWrappers

JLLWrappers.@generate_main_file_header("NetCDF")
JLLWrappers.@generate_main_file("NetCDF", UUID("7243133f-43d8-5620-bbf4-c2c921802cf3"))
end  # module NetCDF_jll

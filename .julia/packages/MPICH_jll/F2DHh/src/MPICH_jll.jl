# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule MPICH_jll
using Base
using Base: UUID
using LazyArtifacts
using MPIPreferences
Base.include(@__MODULE__, joinpath("..", ".pkg", "platform_augmentation.jl"))
import JLLWrappers

JLLWrappers.@generate_main_file_header("MPICH")
JLLWrappers.@generate_main_file("MPICH", UUID("7cb0a576-ebde-5e09-9194-50597f1243b4"))
end  # module MPICH_jll

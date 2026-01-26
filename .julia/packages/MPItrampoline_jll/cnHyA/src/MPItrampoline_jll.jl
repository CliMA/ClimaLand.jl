# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule MPItrampoline_jll
using Base
using Base: UUID
using LazyArtifacts
using MPIPreferences
Base.include(@__MODULE__, joinpath("..", ".pkg", "platform_augmentation.jl"))
import JLLWrappers

JLLWrappers.@generate_main_file_header("MPItrampoline")
JLLWrappers.@generate_main_file("MPItrampoline", UUID("f1f71cc9-e9ae-5b93-9b94-4fe0e1ad3748"))
end  # module MPItrampoline_jll

# Use baremodule to shave off a few KB from the serialized `.ji` file
baremodule x265_jll
using Base
using Base: UUID
import JLLWrappers

JLLWrappers.@generate_main_file_header("x265")
JLLWrappers.@generate_main_file("x265", UUID("dfaa095f-4041-5dcd-9319-2fabd8486b76"))
end  # module x265_jll

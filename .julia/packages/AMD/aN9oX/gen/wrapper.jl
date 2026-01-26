# Script to parse SuiteSparse headers and generate Julia wrappers.
using Clang
using Clang.Generators
using JuliaFormatter

function wrapper(name::String, headers::Vector{String})
  @info "Wrapping $name"

  cd(@__DIR__)
  include_dir = joinpath(pwd(), "include")
  options = load_options(joinpath(@__DIR__, "generator.toml"))
  options["general"]["output_file_path"] = joinpath("..", "src", "wrappers", "$(name).jl")
  options["general"]["library_name"] = "lib" * name

  options["general"]["output_ignorelist"] = [
    "SuiteSparse_config_struct",
    "SuiteSparse_start",
    "SuiteSparse_finish",
    "SuiteSparse_malloc",
    "SuiteSparse_calloc",
    "SuiteSparse_realloc",
    "SuiteSparse_free",
    "SuiteSparse_tic",
    "SuiteSparse_toc",
    "SuiteSparse_time",
    "SuiteSparse_hypot",
    "SuiteSparse_divcomplex",
    "SuiteSparse_version",
    "SuiteSparse_long",
    "SuiteSparse_long_max",
    "SuiteSparse_long_idd",
    "SuiteSparse_long_id",
    "SUITESPARSE_DATE",
    "SUITESPARSE_VERSION",
    "SUITESPARSE_MAIN_VERSION",
    "SUITESPARSE_SUB_VERSION",
    "SUITESPARSE_SUBSUB_VERSION",
    "AMD_DATE",
    "AMD_VERSION",
    "AMD_MAIN_VERSION",
    "AMD_SUB_VERSION",
    "AMD_SUBSUB_VERSION",
    "CAMD_DATE",
    "CAMD_VERSION",
    "CAMD_MAIN_VERSION",
    "CAMD_SUB_VERSION",
    "CAMD_SUBSUB_VERSION",
    "COLAMD_DATE",
    "COLAMD_VERSION",
    "COLAMD_MAIN_VERSION",
    "COLAMD_SUB_VERSION",
    "COLAMD_SUBSUB_VERSION",
    "CCOLAMD_DATE",
    "CCOLAMD_VERSION",
    "CCOLAMD_MAIN_VERSION",
    "CCOLAMD_SUB_VERSION",
    "CCOLAMD_SUBSUB_VERSION",
  ]

  args = get_default_args()
  push!(args, "-I$include_dir")

  ctx = create_context(headers, args, options)
  build!(ctx)

  path = options["general"]["output_file_path"]
  format_file(path, YASStyle())

  text = read(path, String)
  text = replace(text, "Clong" => "SS_Int")
  text = replace(text, "end\n\nconst" => "end\n\n\nconst")
  text = replace(text, "\n\nconst" => "\nconst")
  write(path, text)
  return nothing
end

function main(name::String = "all")
  include = joinpath(pwd(), "include")

  if name == "all" || name == "amd"
    wrapper("amd", [joinpath(include, "amd.h")])
  end

  if name == "all" || name == "camd"
    wrapper("camd", [joinpath(include, "camd.h")])
  end

  if name == "all" || name == "colamd"
    wrapper("colamd", [joinpath(include, "colamd.h")])
  end

  if name == "all" || name == "ccolamd"
    wrapper("ccolamd", [joinpath(include, "ccolamd.h")])
  end
end

# If we want to use the file as a script with `julia wrapper.jl`
if abspath(PROGRAM_FILE) == @__FILE__
  main()
end

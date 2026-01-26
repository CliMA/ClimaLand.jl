using Clang.Generators
using libpng_jll

function rewrite(ex::Expr)
    if Meta.isexpr(ex, :function) && ex.args[1].args[1] == :png_set_write_fn
        # https://github.com/JuliaIO/PNGFiles.jl/pull/27
        ex = quote
            function png_set_write_fn(png_ptr, io_ptr, write_data_fn, output_flush_fn)
                ccall((:png_set_write_fn, libpng), Cvoid, (png_structrp, Any, png_rw_ptr, png_flush_ptr), png_ptr, io_ptr, write_data_fn, output_flush_fn)
            end
        end |> Base.remove_linenums!
        ex.args[]
    else
        ex
    end
end

cd(@__DIR__)

include_dir = joinpath(libpng_jll.artifact_dir, "include") |> normpath

# wrapper generator options
options = load_options(joinpath(@__DIR__, "generator.toml"))

args = get_default_args()
push!(args, "-I$include_dir", "-DPNG_FLOATING_POINT_SUPPORTED")

header_dir = include_dir
headers = [joinpath(header_dir, "png.h")]

# Skip time_t and jmp_buf
@add_def time_t
@add_def jmp_buf

# create context
ctx = create_context(headers, args, options)

# run generator
build!(ctx, BUILDSTAGE_NO_PRINTING)
for node in get_nodes(ctx.dag)
    for (i, expr) in enumerate(node.exprs)
        node.exprs[i] = rewrite(expr)
    end
end
build!(ctx, BUILDSTAGE_PRINTING_ONLY)

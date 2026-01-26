using SafeTestsets

@safetestset "SparseContainers" begin
    include("sparse_containers.jl")
end
@safetestset "Fused incrememnt" begin
    include("fused_increment.jl")
end
@safetestset "Newtons method" begin
    include("test_newtons_method.jl")
end
@safetestset "Single column ARS" begin
    include("single_column_ARS_test.jl")
end
@safetestset "Callbacks" begin
    include("callbacks.jl")
end
@safetestset "Aqua" begin
    include("aqua.jl")
end
@safetestset "Integrator tests" begin
    include("integrator.jl")
end
# These are all run in parallel on buildkite,
# so let's not waste more compilation time.
@safetestset "Algorithm convergence" begin
    include("convergence_lsrk.jl")
    include("tabulate_convergence_orders_multirate.jl")
    include("tabulate_convergence_orders_rosenbrock.jl")
    include("tabulate_convergence_orders_imex_ssp.jl")
    include("tabulate_convergence_orders_imex_ark.jl")
end
@safetestset "Convergence checker unit tests" begin
    include("test_convergence_checker.jl")
end

using SafeTestsets, Test, Pkg
const GROUP = get(ENV, "GROUP", "All")

function activate_downstream_env()
    Pkg.activate("downstream")
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
    Pkg.instantiate()
end

@time begin
    @testset "SciMLOperators" begin
        if GROUP == "All" || GROUP == "Core"
            @time @safetestset "Scalar Operators" begin
                include("scalar.jl")
            end
            @time @safetestset "Basic Operators" begin
                include("basic.jl")
            end
            @time @safetestset "Matrix Operators" begin
                include("matrix.jl")
            end
            @time @safetestset "Function Operator" begin
                include("func.jl")
            end
            @time @safetestset "Full tests" begin
                include("total.jl")
            end
            @time @safetestset "Zygote.jl" begin
                include("zygote.jl")
            end
            @time @safetestset "Copy methods" begin
                include("copy.jl")
            end
        elseif GROUP == "All" || GROUP == "Downstream"
            activate_downstream_env()
            @time @safetestset "AllocCheck" begin
                include("downstream/alloccheck.jl")
            end
        end
    end
end

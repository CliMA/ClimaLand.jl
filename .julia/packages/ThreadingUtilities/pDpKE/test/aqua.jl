import Aqua
@testset "Aqua.jl" begin
    @time Aqua.test_all(ThreadingUtilities)
end

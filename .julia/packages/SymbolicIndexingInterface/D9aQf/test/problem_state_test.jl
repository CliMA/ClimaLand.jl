using SymbolicIndexingInterface
using Test

sys = SymbolCache([:x, :y, :z], [:a, :b, :c], :t)
prob = ProblemState(;
    u = [1.0, 2.0, 3.0], p = [0.1, 0.2, 0.3], t = 0.5, h = Returns(ones(3)))

for (i, sym) in enumerate(variable_symbols(sys))
    @test getsym(sys, sym)(prob) == prob.u[i]
end
for (i, sym) in enumerate(parameter_symbols(sys))
    @test getp(sys, sym)(prob) == prob.p[i]
end
@test getsym(sys, :t)(prob) == prob.t

@test getsym(sys, :(x + a + t))(prob) == 1.6

@test get_history_function(prob) !== nothing

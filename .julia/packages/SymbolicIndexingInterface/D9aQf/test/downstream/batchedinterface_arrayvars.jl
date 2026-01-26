using Symbolics
using SymbolicIndexingInterface

@variables x[1:2] y z

syss = [
    SymbolCache([x..., y]),
    SymbolCache([x[1], y, z])
]
syms = [
    [x, y],
    [x[1], y]
]
probs = [
    ProblemState(; u = [1.0, 2.0, 3.0]),
    ProblemState(; u = [4.0, 5.0, 6.0])
]

bi = BatchedInterface(zip(syss, syms)...)

@test all(isequal.(variable_symbols(bi), [x..., y]))
@test variable_index.((bi,), [x..., y, z]) == [1, 2, 3, nothing]
@test is_variable.((bi,), [x..., y, z]) == Bool[1, 1, 1, 0]
@test associated_systems(bi) == [1, 1, 1]

getter = getsym(bi)
@test (@inferred getter(probs...)) == [1.0, 2.0, 3.0]
buf = zeros(3)
@inferred getter(buf, probs...)
@test buf == [1.0, 2.0, 3.0]

setter! = setsym(bi)
buf .*= 10
setter!(probs..., buf)

@test state_values(probs[1]) == [10.0, 20.0, 30.0]
@test state_values(probs[2]) == [10.0, 30.0, 6.0]

buf ./= 10

setter!(probs[1], 1, buf)
@test state_values(probs[1]) == [1.0, 2.0, 3.0]

@variables a b[1:2] c

syss = [
    SymbolCache([x..., y], [a, b...]),
    SymbolCache([x[1], y, z], [a, b..., c])
]
syms = [
    [x, y, a, b...],
    [x[1], y, b[2], c]
]
probs = [
    ProblemState(; u = [1.0, 2.0, 3.0], p = [0.1, 0.2, 0.3]),
    ProblemState(; u = [4.0, 5.0, 6.0], p = [0.1, 0.4, 0.5, 0.6])
]

bi = BatchedInterface(zip(syss, syms)...)

buf = getsym(bi)(probs...)
buf .*= 100
setter = setsym_oop(bi)
vals = setter(probs..., buf)
@test length(vals) == length(probs)
@test vals[1][1] == [100.0, 200.0, 300.0]
@test vals[1][2] == [10.0, 20.0, 30.0]
@test vals[2][1] == [100.0, 300.0, 6.0]
@test vals[2][2] == [0.1, 0.4, 30.0, 60.0]

buf ./= 10
vals = setter(probs[1], 1, buf)
@test length(vals) == 2
@test vals[1] == [10.0, 20.0, 30.0]
@test vals[2] == [1.0, 2.0, 3.0]

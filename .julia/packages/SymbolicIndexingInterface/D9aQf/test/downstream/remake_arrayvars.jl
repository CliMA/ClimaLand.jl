using ModelingToolkit
using ModelingToolkit: t_nounits as t
using SymbolicIndexingInterface

@variables x(t)[1:2] y(t)
@named sys = ODESystem(Equation[], t, [x, y], [])
sys = complete(sys)

u0 = [1.0, 2.0, 3.0]
newu0 = remake_buffer(sys, u0, [x, y], ([5.0, 6.0], 7.0))
@test newu0 == [5.0, 6.0, 7.0]

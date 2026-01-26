using ModelingToolkit, SymbolicIndexingInterface
using ModelingToolkit: t_nounits as t, D_nounits as D

@variables x(t)[1:2]
@parameters p[1:2, 1:2] q(t)[1:2] r[1:2]

ev = [x[1] ~ 2.0] => [q ~ -ones(2)]
@mtkbuild sys = ODESystem(
    [D(x) ~ p * x + q + r], t, [x], [p, q, r...]; continuous_events = [ev])
@test is_timeseries_parameter(sys, q)
@test !is_timeseries_parameter(sys, p)
@test !is_parameter(sys, r)
@test is_parameter(sys, r[1])
@test is_parameter(sys, r[2])

prob = ODEProblem(
    sys, [x => ones(2)], (0.0, 10.0), [p => ones(2, 2), q => ones(2), r => 2ones(2)])
@test prob.ps[q] ≈ ones(2)
@test prob.ps[p] ≈ ones(2, 2)
@test prob.ps[r] ≈ 2ones(2)
@test prob.ps[p * q] ≈ 2ones(2)

@test getu(sys, p)(prob) ≈ ones(2, 2)
@test getu(sys, r)(prob) ≈ 2ones(2)

prob.ps[p] = 2ones(2, 2)
@test prob.ps[p] ≈ 2ones(2, 2)
prob.ps[q] = 2ones(2)
@test prob.ps[q] ≈ 2ones(2)
prob.ps[r] = ones(2)
@test prob.ps[r] ≈ ones(2)

setter = setp_oop(sys, p)
newp = setter(prob, 3ones(2, 2))
@test getp(sys, p)(newp) ≈ 3ones(2, 2)
setter = setp_oop(sys, r)
newp = setter(prob, 3ones(2))
@test getp(sys, r)(newp) ≈ 3ones(2)

setter = setsym_oop(sys, p)
_, newp = setter(prob, 3ones(2, 2))
@test getp(sys, p)(newp) ≈ 3ones(2, 2)
setter = setsym_oop(sys, r)
_, newp = setter(prob, 3ones(2))
@test getp(sys, r)(newp) ≈ 3ones(2)

@test prob[x] ≈ ones(2)
prob[x] = 2ones(2)
@test prob[x] ≈ 2ones(2)

setu(sys, p)(prob, 4ones(2, 2))
@test prob.ps[p] ≈ 4ones(2, 2)
setu(sys, r)(prob, 4ones(2))
@test prob.ps[r] ≈ 4ones(2)

setter = setsym_oop(sys, x)
newu, newp = setter(prob, 3ones(2))
@test getu(sys, x)(newu) ≈ 3ones(2)

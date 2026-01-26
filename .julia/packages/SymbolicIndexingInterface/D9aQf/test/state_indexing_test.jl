using SymbolicIndexingInterface
using SymbolicIndexingInterface: NotVariableOrParameter

struct FakeIntegrator{S, U, P, T}
    sys::S
    u::U
    p::P
    t::T
end

SymbolicIndexingInterface.symbolic_container(fp::FakeIntegrator) = fp.sys
SymbolicIndexingInterface.state_values(fp::FakeIntegrator) = fp.u
SymbolicIndexingInterface.parameter_values(fp::FakeIntegrator) = fp.p
SymbolicIndexingInterface.current_time(fp::FakeIntegrator) = fp.t

sys = SymbolCache([:x, :y, :z], [:a, :b, :c], [:t])

@test_throws ErrorException getsym(sys, :q)
@test_throws ErrorException setsym(sys, :q)

u = [1.0, 2.0, 3.0]
p = [11.0, 12.0, 13.0]
t = 0.5
fi = FakeIntegrator(sys, copy(u), copy(p), t)
# checking inference for non-concretely typed arrays will always fail
for (sym, val, newval, check_inference) in [(:x, u[1], 4.0, true)
     (:y, u[2], 4.0, true)
     (:z, u[3], 4.0, true)
     (1, u[1], 4.0, true)
     ([:x, :y], u[1:2], 4ones(2), true)
     ([1, 2], u[1:2], 4ones(2), true)
     ((:z, :y), (u[3], u[2]), (4.0, 5.0), true)
     ((3, 2), (u[3], u[2]), (4.0, 5.0), true)
     ([:x, [:y, :z]], [u[1], u[2:3]],
         [4.0, [5.0, 6.0]], false)
     ([:x, 2:3], [u[1], u[2:3]],
         [4.0, [5.0, 6.0]], false)
     ([:x, (:y, :z)], [u[1], (u[2], u[3])],
         [4.0, (5.0, 6.0)], false)
     ([:x, Tuple(2:3)], [u[1], (u[2], u[3])],
         [4.0, (5.0, 6.0)], false)
     ([:x, [:y], (:z,)], [u[1], [u[2]], (u[3],)],
         [4.0, [5.0], (6.0,)], false)
     ([:x, [:y], (3,)], [u[1], [u[2]], (u[3],)],
         [4.0, [5.0], (6.0,)], false)
     ((:x, [:y, :z]), (u[1], u[2:3]),
         (4.0, [5.0, 6.0]), true)
     ((:x, (:y, :z)), (u[1], (u[2], u[3])),
         (4.0, (5.0, 6.0)), true)
     ((1, (:y, :z)), (u[1], (u[2], u[3])),
         (4.0, (5.0, 6.0)), true)
     ((:x, [:y], (:z,)), (u[1], [u[2]], (u[3],)),
         (4.0, [5.0], (6.0,)), true)
     ((a = :x, b = [:x, :y], c = (d = :z, e = :x)),
         (a = u[1], b = u[1:2],
             c = (d = u[3], e = u[1])),
         (a = 4.0, b = [4.0, 5.0],
             c = (d = 6.0, e = 4.0)), true)]
    get = getsym(sys, sym)
    set! = setsym(sys, sym)
    if check_inference
        @inferred get(fi)
    end
    @test get(fi) == val
    if check_inference
        @inferred set!(fi, newval)
    else
        set!(fi, newval)
    end
    @test get(fi) == newval

    new_states = copy(state_values(fi))

    set!(fi, val)
    @test get(fi) == val

    if check_inference
        @inferred get(u)
    end
    @test get(u) == val
    if check_inference
        @inferred set!(u, newval)
    else
        set!(u, newval)
    end
    @test get(u) == newval
    set!(u, val)
    @test get(u) == val

    if sym isa Union{Vector, Tuple} && any(x -> x isa Union{AbstractArray, Tuple}, sym)
        continue
    end

    if !(sym isa NamedTuple)
        setter = setsym_oop(sys, sym)
        svals, pvals = setter(fi, newval)
        @test svals ≈ new_states
        @test pvals == parameter_values(fi)
        @test_throws ArgumentError setter(state_values(fi), newval)
        @test_throws ArgumentError setter(parameter_values(fi), newval)
    end
end

for (sym, val, check_inference) in [
    (:(x + y), u[1] + u[2], true),
    ([:(x + y), :z], [u[1] + u[2], u[3]], false),
    ((:(x + y), :(z + y)), (u[1] + u[2], u[2] + u[3]), false)
]
    get = getsym(sys, sym)
    if check_inference
        @inferred get(fi)
    end
    @test get(fi) == val
end

let fi = fi, sys = sys
    getter = getsym(sys, [])
    @test getter(fi) == []
    getter = getsym(sys, ())
    @test getter(fi) == ()
    sc = SymbolCache(nothing, [:a, :b], :t)
    fi = FakeIntegrator(sys, nothing, [1.0, 2.0], 3.0)
    getter = getsym(sc, [])
    @test getter(fi) == []
    getter = getsym(sc, ())
    @test getter(fi) == ()
end

for (sym, oldval, newval, check_inference) in [(:a, p[1], 4.0, true)
     (:b, p[2], 5.0, true)
     (:c, p[3], 6.0, true)
     ([:a, :b], p[1:2], [4.0, 5.0], true)
     ((:c, :b), (p[3], p[2]), (6.0, 5.0), true)
     ([:x, :a], [u[1], p[1]], [4.0, 5.0], false)
     ((:y, :b), (u[2], p[2]), (5.0, 6.0), true)]
    get = getsym(fi, sym)
    set! = setsym(fi, sym)
    if check_inference
        @inferred get(fi)
    end
    @test get(fi) == oldval
    if check_inference
        @inferred set!(fi, newval)
    else
        set!(fi, newval)
    end
    @test get(fi) == newval

    newu = copy(state_values(fi))
    newp = copy(parameter_values(fi))

    set!(fi, oldval)
    @test get(fi) == oldval

    oop_setter = setsym_oop(sys, sym)
    uvals, pvals = oop_setter(fi, newval)
    @test uvals ≈ newu
    @test pvals ≈ newp
    @test_throws ArgumentError oop_setter(state_values(fi), newval)
    @test_throws ArgumentError oop_setter(parameter_values(fi), newval)
end

for (sym, val, check_inference) in [
    (:t, t, true),
    ([:x, :a, :t], [u[1], p[1], t], false),
    ((:x, :a, :t), (u[1], p[1], t), false)
]
    get = getsym(fi, sym)
    if check_inference
        @inferred get(fi)
    end
    @test get(fi) == val

    @test_throws NotVariableOrParameter setsym_oop(fi, sym)
end

struct FakeSolution{S, U, P, T}
    sys::S
    u::U
    p::P
    t::T
end

SymbolicIndexingInterface.is_timeseries(::Type{<:FakeSolution}) = Timeseries()
function SymbolicIndexingInterface.is_timeseries(::Type{<:FakeSolution{
        S, U, P, Nothing}}) where {S, U, P}
    NotTimeseries()
end
SymbolicIndexingInterface.symbolic_container(fp::FakeSolution) = fp.sys
SymbolicIndexingInterface.state_values(fp::FakeSolution) = fp.u
SymbolicIndexingInterface.parameter_values(fp::FakeSolution) = fp.p
SymbolicIndexingInterface.current_time(fp::FakeSolution) = fp.t

sys = SymbolCache([:x, :y, :z], [:a, :b, :c], [:t])
u = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0], [10.0, 11.0, 12.0]]
t = [1.5, 2.0, 2.3, 4.0]
sol = FakeSolution(sys, u, p, t)

xvals = getindex.(sol.u, 1)
yvals = getindex.(sol.u, 2)
zvals = getindex.(sol.u, 3)

for (sym, ans, check_inference) in [(:x, xvals, true)
     (:y, yvals, true)
     (:z, zvals, true)
     (1, xvals, true)
     ([:x, :y], vcat.(xvals, yvals), true)
     (1:2, vcat.(xvals, yvals), true)
     ([:x, 2], vcat.(xvals, yvals), true)
     ((:z, :y), tuple.(zvals, yvals), true)
     ((3, 2), tuple.(zvals, yvals), true)
     ([:x, [:y, :z]],
         vcat.(xvals, [[x] for x in vcat.(yvals, zvals)]),
         false)
     ([:x, (:y, :z)],
         vcat.(xvals, tuple.(yvals, zvals)), false)
     ([1, (:y, :z)],
         vcat.(xvals, tuple.(yvals, zvals)), false)
     ([:x, [:y, :z], (:x, :z)],
         vcat.(xvals, [[x] for x in vcat.(yvals, zvals)],
             tuple.(xvals, zvals)),
         false)
     ([:x, [:y, 3], (1, :z)],
         vcat.(xvals, [[x] for x in vcat.(yvals, zvals)],
             tuple.(xvals, zvals)),
         false)
     ((:x, [:y, :z]),
         tuple.(xvals, vcat.(yvals, zvals)), true)
     ((:x, (:y, :z)),
         tuple.(xvals, tuple.(yvals, zvals)), true)
     ((:x, [:y, :z], (:z, :y)),
         tuple.(xvals, vcat.(yvals, zvals),
             tuple.(zvals, yvals)),
         true)
     ([:x, :a], vcat.(xvals, p[1]), false)
     ((:y, :b), tuple.(yvals, p[2]), true)
     (:t, t, true)
     ([:x, :a, :t], vcat.(xvals, p[1], t), false)
     ((:x, :a, :t), tuple.(xvals, p[1], t), true)]
    get = getsym(sys, sym)
    if check_inference
        @inferred get(sol)
    end
    @test get(sol) == ans
    for i in [rand(eachindex(u)), CartesianIndex(1), :,
        rand(Bool, length(u)), rand(eachindex(u), 3), 1:3]
        if check_inference
            @inferred get(sol, i)
        end
        @test get(sol, i) == ans[i]
    end
end

for (sym, val, check_inference) in [
    (:(x + y), xvals .+ yvals, true),
    ([:(x + y), :z], vcat.(xvals .+ yvals, zvals), false),
    ((:(x + y), :(z + y)), tuple.(xvals .+ yvals, yvals .+ zvals), false)
]
    get = getsym(sys, sym)
    if check_inference
        @inferred get(sol)
    end
    @test get(sol) == val
    for i in [rand(eachindex(u)), CartesianIndex(1), :,
        rand(Bool, length(u)), rand(eachindex(u), 3), 1:3]
        if check_inference
            @inferred get(sol, i)
        end
        @test get(sol, i) == val[i]
    end
end

for (sym, val) in [(:a, p[1])
                   (:b, p[2])
                   (:c, p[3])
                   ([:a, :b], p[1:2])
                   ((:c, :b), (p[3], p[2]))]
    get = getsym(sys, sym)
    @inferred get(sol)
    @test get(sol) == val
end

let sol = sol, sys = sys
    getter = getsym(sys, [])
    @test getter(sol) == [[] for _ in 1:length(sol.t)]
    getter = getsym(sys, ())
    @test getter(sol) == [() for _ in 1:length(sol.t)]
    sc = SymbolCache(nothing, [:a, :b], :t)
    sol = FakeSolution(sys, [], [1.0, 2.0], [])
    getter = getsym(sc, [])
    @test getter(sol) == []
    getter = getsym(sc, ())
    @test getter(sol) == []
end

sys = SymbolCache([:x, :y, :z], [:a, :b, :c])
u = [1.0, 2.0, 3.0]
p = [10.0, 20.0, 30.0]
fs = FakeSolution(sys, u, p, nothing)
@test is_timeseries(fs) == NotTimeseries()

for (sym, val, check_inference) in [
    (:x, u[1], true),
    (1, u[1], true),
    ([:x, :y], u[1:2], true),
    ((:x, :y), Tuple(u[1:2]), true),
    (1:2, u[1:2], true),
    ([:x, 2], u[1:2], true),
    ((:x, 2), Tuple(u[1:2]), true),
    ([1, 2], u[1:2], true),
    ((1, 2), Tuple(u[1:2]), true),
    (:a, p[1], true),
    ([:a, :b], p[1:2], true),
    ((:a, :b), Tuple(p[1:2]), true),
    ([:x, :a], [u[1], p[1]], false),
    ((:x, :a), (u[1], p[1]), true),
    ([1, :a], [u[1], p[1]], false),
    ((1, :a), (u[1], p[1]), true),
    (:(x + y + a + b), u[1] + u[2] + p[1] + p[2], true),
    ([:(x + a), :(y + b)], [u[1] + p[1], u[2] + p[2]], true),
    ((:(x + a), :(y + b)), (u[1] + p[1], u[2] + p[2]), true)
]
    getter = getsym(sys, sym)
    if check_inference
        @inferred getter(fs)
    end
    @test getter(fs) == val
end

struct NonMarkovianWrapper{S <: SymbolCache}
    sys::S
end

SymbolicIndexingInterface.symbolic_container(hw::NonMarkovianWrapper) = hw.sys
SymbolicIndexingInterface.is_markovian(::NonMarkovianWrapper) = false
function SymbolicIndexingInterface.observed(hw::NonMarkovianWrapper, sym)
    let inner = observed(hw.sys, sym)
        fn(u, h, p, t) = inner(u .+ h(t - 0.1), p, t)
    end
end
function SymbolicIndexingInterface.get_history_function(fs::FakeSolution)
    t -> t .* ones(length(fs.u[1]))
end
function SymbolicIndexingInterface.get_history_function(fi::FakeIntegrator)
    t -> t .* ones(length(fi.u))
end

sys = NonMarkovianWrapper(SymbolCache([:x, :y, :z], [:a, :b, :c], :t))
u0 = [1.0, 2.0, 3.0]
u = [u0 .* i for i in 1:11]
p = [10.0, 20.0, 30.0]
ts = 0.0:0.1:1.0

fi = FakeIntegrator(sys, u0, p, ts[1])
fs = FakeSolution(sys, u, p, ts)
getter = getsym(sys, :(x + y))
@test getter(fi) ≈ 2.8
@test getter(fs) ≈ [3.0i + 2(ts[i] - 0.1) for i in 1:11]
@test getter(fs, 1) ≈ 2.8

pstate = ProblemState(; u = u0, p = p, t = ts[1], h = t -> t .* ones(length(u0)))
@test getter(pstate) ≈ 2.8

struct TupleObservedWrapper{S}
    sys::S
end
SymbolicIndexingInterface.symbolic_container(t::TupleObservedWrapper) = t.sys
SymbolicIndexingInterface.supports_tuple_observed(::TupleObservedWrapper) = true

@testset "Tuple observed" begin
    sc = SymbolCache([:x, :y, :z], [:a, :b, :c])
    sys = TupleObservedWrapper(sc)
    ps = ProblemState(; u = [1.0, 2.0, 3.0], p = [0.1, 0.2, 0.3])
    getter = getsym(sys, (:(x + y), :(y + z)))
    @test all(getter(ps) .≈ (3.0, 5.0))
    @test getter(ps) isa Tuple
    @test_nowarn @inferred getter(ps)
    getter = getsym(sys, (:(a + b), :(b + c)))
    @test all(getter(ps) .≈ (0.3, 0.5))
    @test getter(ps) isa Tuple
    @test_nowarn @inferred getter(ps)

    sc = SymbolCache([:x, :y, :z], [:a, :b, :c], :t)
    sys = TupleObservedWrapper(sc)
    ps = ProblemState(; u = [1.0, 2.0, 3.0], p = [0.1, 0.2, 0.3], t = 0.1)
    getter = getsym(sys, (:(x + y), :(y + t)))
    @test all(getter(ps) .≈ (3.0, 2.1))
    @test getter(ps) isa Tuple
    @test_nowarn @inferred getter(ps)
    getter = getsym(sys, (:(a + b), :(b + c)))
    @test all(getter(ps) .≈ (0.3, 0.5))
    @test getter(ps) isa Tuple
    @test_nowarn @inferred getter(ps)
end

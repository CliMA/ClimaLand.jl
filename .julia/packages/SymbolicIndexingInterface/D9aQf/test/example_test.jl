using SymbolicIndexingInterface
using Test

struct SystemMockup
    static::Bool
    vars::Vector{Symbol}
    params::Vector{Symbol}
    indepvar::Union{Symbol, Nothing}
end

SymbolicIndexingInterface.is_variable(sys::SystemMockup, sym) = sym in sys.vars
function SymbolicIndexingInterface.variable_index(sys::SystemMockup, sym, t = nothing)
    if !constant_structure(sys) && t === nothing
        error("time index must be present")
    end
    findfirst(isequal(sym), variable_symbols(sys, t))
end
function SymbolicIndexingInterface.variable_symbols(sys::SystemMockup, i = nothing)
    return constant_structure(sys) ? sys.vars : circshift(sys.vars, i)
end
SymbolicIndexingInterface.is_parameter(sys::SystemMockup, sym) = sym in sys.params
function SymbolicIndexingInterface.parameter_index(sys::SystemMockup, sym)
    findfirst(isequal(sym), sys.params)
end
SymbolicIndexingInterface.parameter_symbols(sys::SystemMockup) = sys.params
function SymbolicIndexingInterface.is_independent_variable(sys::SystemMockup, sym)
    sys.indepvar !== nothing && isequal(sym, sys.indepvar)
end
function SymbolicIndexingInterface.independent_variable_symbols(sys::SystemMockup)
    if sys.indepvar === nothing
        return []
    else
        return [sys.indepvar]
    end
end
function SymbolicIndexingInterface.is_observed(sys::SystemMockup, sym)
    is_variable(sys, sym) || is_parameter(sys, sym) || is_independent_variable(sys, sym)
end
function SymbolicIndexingInterface.observed(sys::SystemMockup, sym, states = nothing)
    if !constant_structure(sys) && states === nothing
        error("States required")
    end
    states = states isa Vector ? states : variable_symbols(sys, states)
    if is_variable(sys, sym)
        return is_time_dependent(sys) ?
               (u, p, t) -> u[findfirst(isequal(sym), states)] :
               (u, p) -> u[findfirst(isequal(sym), states)]
    end
    idx = parameter_index(sys, sym)
    if idx !== nothing
        return is_time_dependent(sys) ? (u, p, t) -> p[idx] : (u, p) -> p[idx]
    end
    if is_independent_variable(sys, sym)
        return is_time_dependent(sys) ? (u, p, t) -> t : (u, p) -> 1
    end
end
SymbolicIndexingInterface.is_time_dependent(sys::SystemMockup) = isequal(sys.indepvar, :t)
SymbolicIndexingInterface.constant_structure(sys::SystemMockup) = sys.static
SymbolicIndexingInterface.all_variable_symbols(sys::SystemMockup) = sys.vars
function SymbolicIndexingInterface.all_symbols(sys::SystemMockup)
    vcat(sys.vars, sys.params, independent_variable_symbols(sys))
end

sys = SystemMockup(true, [:x, :y, :z], [:a, :b, :c], :t)

@test all(is_variable.((sys,), [:x, :y, :z]))
@test all(.!is_variable.((sys,), [:a, :b, :c, :t, :p, :q, :r]))
@test all(variable_index.((sys,), [:x, :z, :y]) .== [1, 3, 2])
@test all(variable_index.((sys,), [:a, :b, :c, :t, :p, :q, :r]) .=== nothing)
@test all(is_parameter.((sys,), [:a, :b, :c]))
@test all(.!is_parameter.((sys,), [:x, :y, :z, :t, :p, :q, :r]))
@test all(parameter_index.((sys,), [:c, :a, :b]) .== [3, 1, 2])
@test all(parameter_index.((sys,), [:x, :y, :z, :t, :p, :q, :r]) .=== nothing)
@test all(.!is_timeseries_parameter.((sys,), [:x, :y, :z, :t, :p, :q, :r])) # fallback even if not implemented
@test all(timeseries_parameter_index.((sys,), [:x, :y, :z, :t, :p, :q, :r]) .=== nothing) # fallback
@test is_independent_variable(sys, :t)
@test all(.!is_independent_variable.((sys,), [:x, :y, :z, :a, :b, :c, :p, :q, :r]))
@test all(is_observed.((sys,), [:x, :y, :z, :a, :b, :c, :t]))
@test all(observed(sys, :x)(1:3, 4:6, 1.5) .== 1)
@test all(observed(sys, :y)(1:3, 4:6, 1.5) .== 2)
@test all(observed(sys, :z)(1:3, 4:6, 1.5) .== 3)
@test observed(sys, :a)(1:3, 4:6, 1.5) == 4
@test observed(sys, :b)(1:3, 4:6, 1.5) == 5
@test observed(sys, :c)(1:3, 4:6, 1.5) == 6
@test observed(sys, :t)(1:3, 4:6, 1.5) == 1.5
@test is_time_dependent(sys)
@test constant_structure(sys)
@test variable_symbols(sys) == [:x, :y, :z]
@test parameter_symbols(sys) == [:a, :b, :c]
@test independent_variable_symbols(sys) == [:t]
@test all_variable_symbols(sys) == [:x, :y, :z]
@test sort(all_symbols(sys)) == [:a, :b, :c, :t, :x, :y, :z]
@test default_values(sys) == Dict() # fallback even if not implemented

sys = SystemMockup(true, [:x, :y, :z], [:a, :b, :c], nothing)

@test !is_time_dependent(sys)
@test all(observed(sys, :x)(1.0:3.0, 4:6) .== 1.0)
@test all(observed(sys, :y)(1.0:3.0, 4:6) .== 2.0)
@test all(observed(sys, :z)(1.0:3.0, 4:6) .== 3.0)
@test observed(sys, :a)(1:3, 4:6) == 4
@test observed(sys, :b)(1:3, 4:6) == 5
@test observed(sys, :c)(1:3, 4:6) == 6
@test constant_structure(sys)
@test variable_symbols(sys) == [:x, :y, :z]
@test parameter_symbols(sys) == [:a, :b, :c]
@test independent_variable_symbols(sys) == []
@test all_variable_symbols(sys) == [:x, :y, :z]
@test sort(all_symbols(sys)) == [:a, :b, :c, :x, :y, :z]

sys = SystemMockup(false, [:x, :y, :z], [:a, :b, :c], :t)
@test !constant_structure(sys)
for variable in [:x, :y, :z, :a, :b, :c, :t]
    @test_throws ErrorException variable_index(sys, variable)
    @test_throws ErrorException observed(sys, variable)
end
@test all(variable_index.((sys,), [:z, :y, :x], 1) .== [1, 3, 2])
@test all(variable_index.((sys,), [:a, :b, :c, :t], 1) .== nothing)
@test all(observed(sys, :x, 2)(1:3, 4:6, 1.5) .== 3)
@test all(observed(sys, :y, 2)(1:3, 4:6, 1.5) .== 1)
@test all(observed(sys, :z, 2)(1:3, 4:6, 1.5) .== 2)
@test observed(sys, :a, 2)(1:3, 4:6, 1.5) == 4
@test observed(sys, :b, 2)(1:3, 4:6, 1.5) == 5
@test observed(sys, :c, 2)(1:3, 4:6, 1.5) == 6
@test observed(sys, :t, 2)(1:3, 4:6, 1.5) == 1.5
@test_throws Exception variable_symbols(sys)
@test variable_symbols(sys, 1) == [:z, :x, :y]
@test variable_symbols(sys, 2) == [:y, :z, :x]
@test variable_symbols(sys, 3) == [:x, :y, :z]
@test parameter_symbols(sys) == [:a, :b, :c]
@test independent_variable_symbols(sys) == [:t]
@test all_variable_symbols(sys) == [:x, :y, :z]
@test sort(all_symbols(sys)) == [:a, :b, :c, :t, :x, :y, :z]

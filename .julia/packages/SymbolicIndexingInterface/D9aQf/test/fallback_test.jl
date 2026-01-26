using SymbolicIndexingInterface
using Test

struct Wrapper{W}
    wrapped::W
end

SymbolicIndexingInterface.symbolic_container(w::Wrapper) = w.wrapped

sc = SymbolCache([:x, :y, :z], [:a, :b], [:t])
sys = Wrapper(sc)

all_syms = [:x, :y, :z, :a, :b, :t]
@test is_variable.((sys,), all_syms) == is_variable.((sc,), all_syms)
@test variable_index.((sys,), all_syms) == variable_index.((sc,), all_syms)
@test is_parameter.((sys,), all_syms) == is_parameter.((sc,), all_syms)
@test parameter_index.((sys,), all_syms) == parameter_index.((sc,), all_syms)
@test is_independent_variable.((sys,), all_syms) ==
      is_independent_variable.((sc,), all_syms)
@test is_observed.((sys,), all_syms) == is_observed.((sc,), all_syms)
@test is_time_dependent(sys) == is_time_dependent(sc)
@test constant_structure(sys) == constant_structure(sc)
@test variable_symbols(sys) == variable_symbols(sc)
@test parameter_symbols(sys) == parameter_symbols(sc)
@test independent_variable_symbols(sys) == independent_variable_symbols(sc)
@test all_variable_symbols(sys) == variable_symbols(sc)
@test all_symbols(sys) == all_symbols(sc)

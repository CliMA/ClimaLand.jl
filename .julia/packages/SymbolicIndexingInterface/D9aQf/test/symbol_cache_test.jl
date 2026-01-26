using SymbolicIndexingInterface
using Test

sc = SymbolCache(
    [:x, :y, :z], [:a, :b], [:t]; defaults = Dict(:x => 1, :y => :(2b), :b => :(2a + x)))

@test all(is_variable.((sc,), [:x, :y, :z]))
@test all(.!is_variable.((sc,), [:a, :b, :t, :q]))
@test variable_index.((sc,), [:x, :y, :z, :a]) == [1, 2, 3, nothing]
@test all(is_parameter.((sc,), [:a, :b]))
@test all(.!is_parameter.((sc,), [:x, :y, :z, :t, :q]))
@test parameter_index.((sc,), [:a, :b, :x]) == [1, 2, nothing]
@test is_independent_variable(sc, :t)
@test all(.!is_independent_variable.((sc,), [:x, :y, :z, :a, :b, :q]))
@test all(.!is_observed.((sc,), [:x, :y, :z, :a, :b, :t, :q]))
@test is_time_dependent(sc)
@test constant_structure(sc)
@test variable_symbols(sc) == [:x, :y, :z]
@test sort(parameter_symbols(sc)) == [:a, :b]
@test independent_variable_symbols(sc) == [:t]
@test all_variable_symbols(sc) == [:x, :y, :z]
@test sort(sort(all_symbols(sc))) == [:a, :b, :t, :x, :y, :z]
@test default_values(sc)[:x] == 1
@test default_values(sc)[:y] == :(2b)
@test default_values(sc)[:b] == :(2a + x)

@test symbolic_evaluate(:x, default_values(sc)) == 1
@test symbolic_evaluate(:y, default_values(sc)) == :(2 * (2a + 1))
@test symbolic_evaluate(:(x + y), merge(default_values(sc), Dict(:a => 2))) == 11

@test is_observed(sc, :(x + a + t))
obsfn = observed(sc, :(x + a + t))
@test obsfn(ones(3), 2ones(2), 3.0) == 6.0
obsfn2 = observed(sc, :(x + a + t))
@test obsfn === obsfn2

@test is_observed(sc, [:(x + a), :(a + t)])
obsfn3 = observed(sc, [:(x + a), :(a + t)])
@test obsfn3(ones(3), 2ones(2), 3.0) ≈ [3.0, 5.0]
@test is_observed(sc, [:(x + a) :(y + b); :(x + y) :(a + b)])
obsfn4 = observed(sc, [:(x + a) :(y + b); :(x + y) :(a + b)])
@test size(obsfn4(ones(3), 2ones(2), 3.0)) == (2, 2)
@test obsfn4(ones(3), 2ones(2), 3.0) ≈ [3.0 3.0; 2.0 4.0]
@test is_observed(sc, (:(x + a), :(y + b)))
obsfn5 = observed(sc, (:(x + a), :(y + b)))
@test all(obsfn5(ones(3), 2ones(2), 3.0) .≈ (3.0, 3.0))

@test_throws TypeError observed(sc, [:(x + a), 2])
@test_throws TypeError observed(sc, (:(x + a), 2))

pobsfn1 = parameter_observed(sc, :(a + b + t))
@test pobsfn1(2ones(2), 3.0) == 7.0
pobsfn2 = parameter_observed(sc, [:(a + b + t), :(a + t)])
@test pobsfn2(2ones(2), 3.0) == [7.0, 5.0]
buffer = zeros(2)
pobsfn2(buffer, 2ones(2), 3.0)
@test buffer == [7.0, 5.0]
pobsfn3 = parameter_observed(sc, (:(a + b + t), :(a + t)))
@test pobsfn3(2ones(2), 3.0) == (7.0, 5.0)
buffer = zeros(2)
pobsfn3(buffer, 2ones(2), 3.0)
@test buffer == [7.0, 5.0]

@test_throws TypeError parameter_observed(sc, [:(a + b), 4])
@test_throws TypeError parameter_observed(sc, (:(a + b), 4))

sc = SymbolCache([:x, :y], [:a, :b, :c], :t;
    timeseries_parameters = Dict(
        :b => ParameterTimeseriesIndex(1, 1), :c => ParameterTimeseriesIndex(2, 1)))
@test only(get_all_timeseries_indexes(sc, :(a + c))) == 2
@test only(get_all_timeseries_indexes(sc, [:a, :c])) == 2
@test isempty(get_all_timeseries_indexes(sc, :(2a)))
@test isempty(get_all_timeseries_indexes(sc, [:(2a), :(3a)]))
@test sort(collect(get_all_timeseries_indexes(sc, [:b, :c]))) == [1, 2]

@test_throws ArgumentError SymbolCache([:x, :y], [:a, :b], :t;
    timeseries_parameters = Dict(:c => ParameterTimeseriesIndex(1, 1)))
@test_throws TypeError SymbolCache(
    [:x, :y], [:a, :c], :t; timeseries_parameters = Dict(:c => (1, 1)))
@test_nowarn SymbolCache([:x, :y], [:a, :c], :t;
    timeseries_parameters = Dict(:c => ParameterTimeseriesIndex(1, 1)))

sc = SymbolCache([:x, :y], [:a, :b])
@test !is_time_dependent(sc)
@test sort(all_symbols(sc)) == [:a, :b, :x, :y]
@test is_observed(sc, :(x + b))
obsfn = observed(sc, :(x + b))
@test obsfn(ones(2), 2ones(2)) == 3.0
# make sure the constructor works
@test_nowarn SymbolCache([:x, :y])

@test_throws ArgumentError SymbolCache(
    [:x, :y], [:a, :b]; timeseries_parameters = Dict(:b => ParameterTimeseriesIndex(1, 1)))

sc = SymbolCache()
@test all(.!is_variable.((sc,), [:x, :y, :a, :b, :t]))
@test all(variable_index.((sc,), [:x, :y, :a, :b, :t]) .== nothing)
@test variable_symbols(sc) == []
@test all(.!is_parameter.((sc,), [:x, :y, :a, :b, :t]))
@test all(parameter_index.((sc,), [:x, :y, :a, :b, :t]) .== nothing)
@test parameter_symbols(sc) == []
@test all(.!is_independent_variable.((sc,), [:x, :y, :a, :b, :t]))
@test independent_variable_symbols(sc) == []
@test !is_time_dependent(sc)
@test all_variable_symbols(sc) == []
@test all_symbols(sc) == []
@test isempty(default_values(sc))

sc = SymbolCache(nothing, nothing, :t)
@test all(.!is_independent_variable.((sc,), [:x, :y, :a, :b]))
@test is_independent_variable(sc, :t)
@test independent_variable_symbols(sc) == [:t]
@test is_time_dependent(sc)
@test all_variable_symbols(sc) == []
@test all_symbols(sc) == [:t]
@test isempty(default_values(sc))

sc = SymbolCache(nothing, nothing, [:t1, :t2, :t3])
@test all(is_independent_variable.((sc,), [:t1, :t2, :t3]))
@test independent_variable_symbols(sc) == [:t1, :t2, :t3]

sc2 = copy(sc)
@test sc.variables == sc2.variables
@test sc.parameters == sc2.parameters
@test sc.independent_variables == sc2.independent_variables

sc = SymbolCache()
for sym in [1, :a, :(a + b), "foo", [:a, :b], [:(a + b), :c]]
    @test only(get_all_timeseries_indexes(sc, sym)) == ContinuousTimeseries()
end

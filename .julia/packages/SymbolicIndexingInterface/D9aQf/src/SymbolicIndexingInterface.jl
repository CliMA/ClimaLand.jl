module SymbolicIndexingInterface

using RuntimeGeneratedFunctions
import StaticArraysCore: MArray, similar_type
import ArrayInterface
using Accessors: @reset

RuntimeGeneratedFunctions.init(@__MODULE__)

export ScalarSymbolic, ArraySymbolic, NotSymbolic, symbolic_type, hasname, getname,
       Timeseries, NotTimeseries, is_timeseries, is_parameter_timeseries
include("trait.jl")

export is_variable, variable_index, variable_symbols, is_parameter, parameter_index,
       is_timeseries_parameter, timeseries_parameter_index, ParameterTimeseriesIndex,
       parameter_symbols, is_independent_variable, independent_variable_symbols,
       is_observed, observed, parameter_observed,
       ContinuousTimeseries, get_all_timeseries_indexes,
       is_time_dependent, is_markovian, constant_structure, symbolic_container,
       all_variable_symbols, all_symbols, solvedvariables, allvariables, default_values,
       symbolic_evaluate
include("index_provider_interface.jl")

export SymbolCache
include("symbol_cache.jl")

export parameter_values, set_parameter!, finalize_parameters_hook!,
       get_parameter_timeseries_collection, with_updated_parameter_timeseries_values,
       state_values, set_state!, current_time, get_history_function
include("value_provider_interface.jl")

export ParameterTimeseriesCollection
include("parameter_timeseries_collection.jl")

export getp, setp, setp_oop
include("parameter_indexing.jl")

export getsym, setsym, getu, setu
include("state_indexing.jl")

export BatchedInterface, setsym_oop, associated_systems
include("batched_interface.jl")

export ProblemState
include("problem_state.jl")

export ParameterIndexingProxy, show_params
include("parameter_indexing_proxy.jl")

export remake_buffer
include("remake.jl")

include("despecialize.jl")

end

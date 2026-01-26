"""
    struct ProblemState
    function ProblemState(; u = nothing, p = nothing, t = nothing, h = nothing)

A value provider struct which can be used as an argument to the function returned by
[`getsym`](@ref) or [`setsym`](@ref). It stores the state vector, parameter object and
current time, and forwards calls to [`state_values`](@ref), [`parameter_values`](@ref),
[`current_time`](@ref), [`set_state!`](@ref), [`set_parameter!`](@ref) to the contained
objects.

A history function may be provided using the `h` keyword, which will be returned with
[`get_history_function`](@ref).
"""
struct ProblemState{U, P, T, H}
    u::U
    p::P
    t::T
    h::H
end

function ProblemState(; u = nothing, p = nothing, t = nothing, h = nothing)
    ProblemState(u, p, t, h)
end

state_values(ps::ProblemState) = ps.u
parameter_values(ps::ProblemState) = ps.p
current_time(ps::ProblemState) = ps.t
set_state!(ps::ProblemState, val, idx) = set_state!(ps.u, val, idx)
set_parameter!(ps::ProblemState, val, idx) = set_parameter!(ps.p, val, idx)
get_history_function(ps::ProblemState) = ps.h

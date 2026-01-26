function _string_to_bool(s::AbstractString)
    b = tryparse(Bool, s)
    if b isa Nothing
        return false
    else
        return b::Bool
    end
end

function _is_suppress_warning()
    s = get(ENV, "SUPPRESS_THREADINGUTILITIES_WARNING", "")
    b = _string_to_bool(s)::Bool
    return b
end

function _is_julia_exclusive()
    s = strip(get(ENV, "JULIA_EXCLUSIVE", "0"))
    # From https://github.com/JuliaLang/julia/blob/efc986003356be5daea6ccd9e38008de961c18aa/doc/src/manual/environment-variables.md:
    # If set to anything besides 0, then Julia's thread policy is consistent
    # with running on a dedicated machine: the master thread is on proc 0, and
    # threads are affinitized. Otherwise, Julia lets the operating system handle
    # thread policy.
    i = tryparse(Int, s)
    if isempty(s)
        return false
    elseif s == "0" || i == 0
        return false
    else
        return true
    end
end

function _print_exclusivity_warning()
    is_julia_exclusive = _is_julia_exclusive()::Bool
    if !is_julia_exclusive
        suppress_warning = _is_suppress_warning()::Bool
        if !suppress_warning
            msg = string(
                "The JULIA_EXCLUSIVE environment variable is not set to 1. ",
                "Therefore, the kernel is allowed to move threads to different ",
                "cores. We recommend that you set the JULIA_EXCLUSIVE ",
                "environment variable to 1. To suppress this warning, set the ",
                "SUPPRESS_THREADINGUTILITIES_WARNING environment ",
                "variable to 1.",
            )
            @debug(msg)
        end
    end
    return nothing
end


# default is disable = true, so explicitly enable first
using Preferences, UUIDs
Preferences.set_preferences!(UUID("781d530d-4396-4725-bb49-402e4bee1e77"), "disable" => false)

using Test, TruncatedStacktraces

@testset "Test that VERBOSE can remove the notice message" begin
    TruncatedStacktraces.VERBOSE[] = false
    error_msg = Ref{String}()
    try
        x
    catch e
        io = IOBuffer()
        showerror(io, e)
        error_msg[] = String(take!(io))
    end
    @static if VERSION >= v"1.9.0-rc1"
        actual_error_msg = "UndefVarError: `x` not defined" *
                           "\n\nSome of the types have been truncated in the" *
                           " stacktrace for improved reading. To emit complete " *
                           "information\nin the stack trace, evaluate " *
                           "`TruncatedStacktraces.VERBOSE[] = true` and re-run the code.\n"
    else
        actual_error_msg = "UndefVarError: x not defined" *
                           "\n\nSome of the types have been truncated in the" *
                           " stacktrace for improved reading. To emit complete " *
                           "information\nin the stack trace, evaluate " *
                           "`TruncatedStacktraces.VERBOSE[] = true` and re-run the code.\n"
    end
    # Printing the hint message is broken in Julia 1.6
    @static if v"1.6" <= VERSION < v"1.7"
        @test_broken error_msg[] == actual_error_msg
    else
        @test error_msg[] == actual_error_msg
    end

    TruncatedStacktraces.VERBOSE[] = true
    try
        x
    catch e
        io = IOBuffer()
        showerror(io, e)
        error_msg[] = String(take!(io))
    end
    @static if VERSION >= v"1.9.0-rc1"
        actual_error_msg = "UndefVarError: `x` not defined"
    else
        actual_error_msg = "UndefVarError: x not defined"
    end
    @test error_msg[] == actual_error_msg
end

module TerminalTools

using Dates
import REPL: Terminals
import Base: TTY

struct TimeoutException <: Exception
    timeout::Float64 # seconds
end

"""
    timeout_call(f, timeout; pollint=0.1)
Execute function `f()` with a maximum timeout `timeout` seconds.

An `TimeoutException(timeout)` exception will be thrown if it exceeds
the maximum timeout. If `f()` exits with error, it will be rethrown.
"""
function timeout_call(f::Function, timeout::Real; pollint=0.1)
    start = now()

    t = @task f()
    schedule(t)

    while !istaskdone(t)
        if (now()-start).value >= 1000timeout
            schedule(t, TimeoutException(timeout), error=true)
            sleep(pollint) # wait a while for the task to update its state
            break
        end
        sleep(pollint)
    end

    if t.state == :failed
        throw(t.exception)
    else
        return t.result
    end
end

function with_raw(f, tty::Terminals.TTYTerminal)
    Terminals.raw!(tty, true)
    try
        return f()
    finally
        Terminals.raw!(tty, false)
    end
end

query_terminal(msg, io::IO; kwargs...) = ""
query_terminal(msg, regex, io::IO; kwargs...) = ("", )
function query_terminal(msg, tty::TTY; timeout=1)
    term = Terminals.TTYTerminal("", stdin, tty, stderr)
    try
        timeout_call(timeout; pollint=timeout/100) do
            with_raw(term) do
                write(tty, msg)
                return @static if Sys.iswindows()
                    response = ""
                    while !endswith(response, 'c')
                        response *= read(stdin, Char)
                    end
                    response
                else
                    transcode(String, readavailable(tty))
                end
            end
        end
    catch e
        e isa TimeoutException && return ""
        rethrow()
    end
end
function query_terminal(msg, regex, tty::TTY; kwargs...)
    response = query_terminal(msg, tty; kwargs...)
    m = match(regex, response)
    isnothing(m) ? ("", ) : Tuple(m.captures)
end

# get_window_title(tty=stdout) = query_terminal("\033[21t", r"\033\]l(.*)\033\\\\")
# get_window_icon_title(tty=stdout) = query_terminal("\033[20t", r"\033\]L(.*)\033\\\\")
# get_window_size(tty=stdout) = parse.(Int, query_terminal("\033[14t", r"\033\[4;(?<row>[0-9]*);(?<col>[0-9]*)t"))
# get_window_position(tty=stdout) = parse.(Int, query_terminal("\033[13t", r"\033\[3;(?<row>[0-9]*);(?<col>[0-9]*)t"))
# get_screen_size(tty=stdout) = parse.(Int, query_terminal("\033[19t", r"\033\[9;(?<row>[0-9]*);(?<col>[0-9]*)t"))
# get_text_area(tty=stdout) = parse.(Int, query_terminal("\033[18t", r"\033\[8;(?<row>[0-9]*);(?<col>[0-9]*)t"))
# get_text_area(tty=stdout) = displaysize(tty)
# get_cursor_position(tty=stdout) = parse.(Int, query_terminal("\033[6n", r"\033\[(?<row>[0-9]*);(?<col>[0-9]*)R"))

end # moudle

# Logging

## Overview

Logging is a crucial tool for debugging, monitoring, and understanding the behavior of your applications. ClimaComms extends Julia's built-in logging functionality (provided by [`Logging.jl`](https://docs.julialang.org/en/v1/stdlib/Logging/)) to work seamlessly in distributed computing environments, particularly with MPI.

### Julia's Logging System

Julia's standard library provides a flexible logging system through `Logging.jl`. Key concepts include:

- **Logger**: An object that handles log messages, determining how they are formatted and where they are sent
- **Log Level**: Indicates the severity/importance of a message
- **Log Record**: Contains the message and associated metadata (level, module, line number, etc.)

The default logger in Julia (as of v1.11) is `Logging.ConsoleLogger()`, which prints messages to the console. You can check your current logger using:

```julia
using Logging; current_logger()
```

## ClimaComms Logging Features

ClimaComms builds upon Julia's logging system by providing specialized loggers for distributed computing.

To set any of the loggers below as the global logger, use `Logging.global_logger(logger)`:

```julia
using ClimaComms, Logging

ctx = ClimaComms.context()
logger = ClimaComms.OnlyRootLogger(ctx)
global_logger(logger)
```

### OnlyRootLogger

`OnlyRootLogger(context)` returns a logger that silences non-root processes.
If using MPI, this logger is enabled on the first [`ClimaComms.init`](@ref) call.

### FileLogger

`FileLogger(context, log_dir)` writes to `stdout` and to a file `log_dir/output.log` simultaneously.
If using MPI, the `FileLogger` separates logs by MPI process into different files, so that process `i` will write to the `rank_$i.log` file in the `$log_dir` directory, where `$log_dir` is of your choosing. 
In this case, `$log_dir/output.log` is a [symbolic link](https://en.wikipedia.org/wiki/Symbolic_link) to the root process logs. 
In other words, you can always look at `$log_dir/output.log` for the output. Logging to `stdout` can be disabled by setting the keyword argument `log_stdout = false`.
```julia
using ClimaComms, Logging
ctx = ClimaComms.context()
logger = ClimaComms.FileLogger(ctx, "logs")

with_logger(logger) do
    @warn "Memory usage high"  # Written to rank-specific log file
end
```
This will output the following in both the REPL and `logs/rank_1.log`:
```julia
┌ Warning: Memory usage high
└ @ Main REPL[6]:2
```


### MPILogger

`MPILogger(context)` adds an MPI rank prefix to all log messages:

```julia
using ClimaComms, Logging
ctx = ClimaComms.context()
logger = MPILogger(ctx)

with_logger(logger) do
    @info "Processing data..."  # Output: [P1]  Info: Processing data...
end
```

## Log Levels

Julia provides four standard log levels, in order of increasing severity:

1. `Debug`: Detailed information for debugging
2. `Info`: General information about program execution
3. `Warn`: Warnings about potential issues
4. `Error`: Error conditions that might still allow the program to continue

See [Julia documentation](https://docs.julialang.org/en/v1/stdlib/Logging/#Log-event-structure) for more detailed information.

You can define custom log levels using `LogLevel`:

```julia
const Trace = LogLevel(-1000)  # Lower number = less severe
@macroexpand @logmsg Trace "Very detailed trace message"
```

To disable all log messages at log levels equal to or less than a given `LogLevel`, use [`Logging.disable_logging(level)`](https://docs.julialang.org/en/v1/stdlib/Logging/#Logging.disable_logging).

## Filtering Log Messages

[LoggingExtras.jl](https://github.com/JuliaLogging/LoggingExtras.jl) provides powerful filtering capabilities through the `EarlyFilteredLogger(filter, logger)`, which takes two arguments:

- `filter(log_args)` is a function which takes in `log_args` and returns a Boolean determining if the message should be logged. `log_args` is a NamedTuple with fields `level`, `_module`, `id` and `group`. Example `filter` functions are provided below in the Common Use Cases.
- `logger` is any existing logger, such as `Logging.ConsoleLogger()` or `MPILogger(ctx)`.

### Common Use Cases

#### How do I save log output to stdout and a file simultaneously?

`ClimaComms.FileLogger` logs to files and stdout simultaneously.
For full customization, `LoggingExtras.TeeLogger(loggers)` composes multiple loggers, allowing for multiple loggers at once as shown below. 

```julia
using Logging, LoggingExtras

io = open("simulation.log", "w")
loggers = (Logging.ConsoleLogger(stdout), Logging.ConsoleLogger(io))
tee_logger = LoggingExtras.TeeLogger(loggers)

with_logger(tee_logger) do
    @warn "Log to stdout and file"
end

close(io)
```

#### How do I filter out warning messages?

```@example
using Logging, LoggingExtras

function no_warnings(log_args)
    return log_args.level != Logging.Warn
end

filtered_logger = EarlyFilteredLogger(no_warnings, Logging.current_logger())

with_logger(filtered_logger) do
    @warn "Hide this warning"
    @info "Display this message"
end
```

#### How do I filter out messages from certain modules?

We can create a custom filter that returns `false` if a log message originates from a list of excluded modules.

The same pattern can be reversed to filter messages only coming from certain modules.
```julia
using Logging, LoggingExtras

module_filter(excluded_modules) = log_args ->
    !(log_args._module in excluded_modules)

ModuleFilteredLogger(excluded) =
    EarlyFilteredLogger(module_filter(excluded), Logging.current_logger())
# To test this logger:
module TestModule
    using Logging
    function log_something()
        @info "This message will appear"
    end
end

excluded = (Main, Base)
with_logger(ModuleFilteredLogger(excluded)) do
    @info "Hide this message"
    TestModule.log_something()
end
```
```julia
[ Info: This message will appear
```

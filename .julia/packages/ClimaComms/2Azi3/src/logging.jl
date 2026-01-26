import Logging, LoggingExtras

export MPILogger, MPIFileLogger

"""
    OnlyRootLogger()
    OnlyRootLogger(ctx::AbstractCommsContext)

Return a logger that silences non-root processes.
    
If no context is passed, obtain the default context via [`context`](@ref).
"""
OnlyRootLogger() = OnlyRootLogger(context())

function OnlyRootLogger(ctx::AbstractCommsContext)
    if iamroot(ctx)
        return Logging.ConsoleLogger()
    else
        return Logging.NullLogger()
    end
end

"""
    MPILogger(context::AbstractCommsContext)
    MPILogger(iostream, context)
    
Add a rank prefix before log messages.

Outputs to `stdout` if no IOStream is given.
"""
MPILogger(ctx::AbstractCommsContext) = MPILogger(stdout, ctx)

function MPILogger(iostream, ctx::AbstractCommsContext)
    pid = mypid(ctx)

    function format_log(io, log)
        print(io, "[P$pid] ")
        println(io, " $(log.level): $(log.message)")
    end

    return LoggingExtras.FormatLogger(format_log, iostream)
end

"""
    FileLogger(context, log_dir; log_stdout = true, min_level = Logging.Info)

Log MPI ranks to different files within the `log_dir`.

The minimum logging level is set using `min_level`.
If `log_stdout = true`, root process logs will be sent to stdout as well.
"""
function FileLogger(
    ctx::AbstractCommsContext,
    log_dir;
    log_stdout = true,
    min_level::Logging.LogLevel = Logging.Info,
)
    return FileLogger(stdout, ctx, log_dir; log_stdout, min_level)
end

function FileLogger(
    io::IO,
    ctx::MPICommsContext,
    log_dir::AbstractString;
    log_stdout = true,
    min_level::Logging.LogLevel = Logging.Info,
)
    mpi_log_dir = joinpath(log_dir, "logs")
    !isdir(mpi_log_dir) && mkpath(mpi_log_dir)
    ClimaComms.barrier(ctx)  # Ensure that the folder is created
    rank = mypid(ctx)
    filepath = abspath(joinpath(mpi_log_dir, "rank_$rank.log"))

    state = Dict(:logger_used => false)

    function min_level_filter(log_args)
        if iamroot(ctx) && !state[:logger_used]
            state[:logger_used] = true
            symlink_path = abspath(joinpath(log_dir, "output.log"))
            symlink(filepath, symlink_path)
        end
        return log_args.level >= min_level
    end

    file_logger = LoggingExtras.FormatLogger(
        format_log,
        filepath,
        append = true,
        always_flush = true,
    )

    filtered_logger =
        LoggingExtras.EarlyFilteredLogger(min_level_filter, file_logger)

    if iamroot(ctx) && log_stdout
        return LoggingExtras.TeeLogger((
            Logging.ConsoleLogger(io),
            filtered_logger,
        ))
    else
        return filtered_logger
    end
end

function FileLogger(
    io::IO,
    ctx::SingletonCommsContext,
    log_dir::AbstractString;
    log_stdout = true,
    min_level::Logging.LogLevel = Logging.Info,
)
    !isdir(log_dir) && mkpath(log_dir)
    filepath = joinpath(log_dir, "output.log")

    function min_level_filter(log_args)
        return log_args.level >= min_level
    end

    file_logger = LoggingExtras.FormatLogger(
        format_log,
        filepath,
        append = true,
        always_flush = true,
    )

    filtered_logger =
        LoggingExtras.EarlyFilteredLogger(min_level_filter, file_logger)

    if log_stdout
        return LoggingExtras.TeeLogger((
            Logging.ConsoleLogger(io),
            filtered_logger,
        ))
    else
        return filtered_logger
    end
end

"""
    format_log(io, args)

Format log messages similarly to `Logging.ConsoleLogger` for the FileLogger

Add box decorations for multiline strings, indentation, and bolding.
"""
function format_log(io::IO, args)
    msg = string(args.message)
    msglines = split(msg, '\n')
    if !isempty(args.kwargs)
        for (key, val) in args.kwargs
            push!(msglines, "$key = $val")
        end
    end

    level_prefix = string(args.level)

    for (i, msg) in enumerate(msglines)
        # Set up the box decoration for multi-line strings
        boxstr =
            length(msglines) == 1 ? "[ " :
            i == 1 ? "┌ " : i < length(msglines) ? "│ " : "└ "
        printstyled(io, boxstr, bold = true)
        if i == 1
            printstyled(io, level_prefix, ": ", bold = true)
        end
        indent = i == 1 ? 0 : 2
        print(io, " "^indent, msg)
        println(io)
    end
end

"""
    with_tempdir(f::Function)

Call `f` on a temporary directory.
"""
function with_tempdir(f::Function, ctx)
    temp_dir = ClimaComms.iamroot(ctx) ? mktempdir() : nothing
    temp_dir = ClimaComms.bcast(ctx, temp_dir)
    return f(temp_dir)
end

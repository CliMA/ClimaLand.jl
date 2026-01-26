using Logging: AbstractLogger, Logging

struct ScopePayloadLogger <: AbstractLogger
    logger::AbstractLogger
    scope::Union{Scope, Nothing}
end

function current_scope()
    logger = Logging.current_logger()
    if logger isa ScopePayloadLogger
        return logger.scope
    end
    return nothing
end

function enter_scope(f, scope)
    return Logging.with_logger(f, ScopePayloadLogger(current_logger(), scope))
end

# Forward actual logging interface:
Logging.handle_message(payload::ScopePayloadLogger, args...; kwargs...) =
    Logging.handle_message(payload.logger, args...; kwargs...)
Logging.shouldlog(payload::ScopePayloadLogger, args...) =
    Logging.shouldlog(payload.logger, args...)
Logging.min_enabled_level(payload::ScopePayloadLogger, args...) =
    Logging.min_enabled_level(payload.logger, args...)
Logging.catch_exceptions(payload::ScopePayloadLogger, args...) =
    Logging.catch_exceptions(payload.logger, args...)

"""
    ScopedValues.with_logger(f, logger::AbstractLogger)

Like `Logging.with_logger` but properly propagate the scope.
"""
function with_logger(f, logger::AbstractLogger)
    @nospecialize
    cpl = Logging.current_logger()
    if cpl isa ScopePayloadLogger
        scope = cpl.scope
    else
        scope = Scope(nothing)
    end
    return Logging.with_logger(f, ScopePayloadLogger(logger, scope))
end

"""
    ScopedValues.current_logger() -> logger::AbstractLogger

Like `Logging.current_logger` but unwraps `ScopePayloadLogger`.
"""
function current_logger()
    logger = Logging.current_logger()
    if logger isa ScopePayloadLogger
        return logger.logger
    else
        return logger
    end
end

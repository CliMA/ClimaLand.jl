"""
    AbstractCFDateTime{T,Torigintuple}

Supertype for all DateTime structures following the CF conventions
where the time instance is represented as a duration of type `T` since a time
origin specified by the value type of the tuple `Torigintuple`.

The tuple is composed of the integers representing the year,
month and day as well as smaller time divisions if necessary.

The type parameter `T` and `Torigintuple` are considered as internal API.
Only reducing the number of type parameters of `AbstractCFDateTime` will be
considered as a breaking change.
"""
abstract type AbstractCFDateTime{T, Torigintuple} <: Dates.TimeType
end


"""
    Period{T,Tfactor,Texponent}

Period wraps a number duration of type `T` where

    duration * factor * 10^exponent

represents the time in seconds
"""
struct Period{T, Tfactor, Texponent}
    duration::T
end


for (CFDateTime, calendar) in [
        (:DateTimeStandard, "standard"),
        (:DateTimeJulian, "julian"),
        (:DateTimeProlepticGregorian, "prolepticgregorian"),
        (:DateTimeAllLeap, "allleap"),
        (:DateTimeNoLeap, "noleap"),
        (:DateTime360Day, "360day"),
        (:DateTimeUTC, "utc"),
        (:DateTimeTAI, "tai"),

    ]
    @eval begin
        struct $CFDateTime{T, Torigintuple} <: AbstractCFDateTime{T, Torigintuple}
            instant::T
        end

    end
end

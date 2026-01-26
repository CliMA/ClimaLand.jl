# year upto attosecond

for (i, func) in enumerate(TIME_NAMES)
    modname = (i <= 9 ? "Dates" : "CFTime")
    @eval begin
        """
            $($modname).$($func)(dt::AbstractCFDateTime) -> Int64

        Extract the $($func) part of an `AbstractCFDateTime` as an `Int64`.
        """
        function $func(dt::AbstractCFDateTime)
            t = _datetuple(dt)
            if length(t) >= $i
                return Int64(t[$i])
            else
                return Int64(0)
            end
        end
    end
end

# Day, Hour, ... Nanosecond constructors from Dates.
# There is no Picosecond constructor in Dates.

for (i, name) in enumerate(TIME_NAMES[1:9])
    function_name = Symbol(uppercasefirst(String(name)))

    @eval begin
        """
            Dates.$($name)(dt::AbstractCFDateTime) -> $($function_name)

        The $($name) part of an `AbstractCFDateTime` as an `$($function_name)`.
        """
        @inline function $function_name(dt::AbstractCFDateTime)
            return $function_name(Dates.$name(dt)) # years and months are special
        end
    end
end

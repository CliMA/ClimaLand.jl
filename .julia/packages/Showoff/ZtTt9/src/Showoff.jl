module Showoff

using Dates

if isdefined(Base, :Ryu)
    include("ryu.jl")
else
    include("grisu.jl")
end

export showoff

# suppress compile errors when there isn't a grisu_ccall macro
macro grisu_ccall(x, mode, ndigits)
    quote end
end

# Fallback
function showoff(xs::AbstractArray, style=:none)
    result = Vector{String}(undef, length(xs))
    buf = IOBuffer()
    for (i, x) in enumerate(xs)
        show(buf, x)
        result[i] = String(take!(buf))
    end

    return result
end


# Floating-point

function concrete_minimum(xs)
    if isempty(xs)
        throw(ArgumentError("argument must not be empty"))
    end

    x_min = first(xs)
    for x in xs
        if isa(x, AbstractFloat) && isfinite(x)
            x_min = x
            break
        end
    end

    for x in xs
        if isa(x, AbstractFloat) && isfinite(x) && x < x_min
            x_min = x
        end
    end
    return x_min
end


function concrete_maximum(xs)
    if isempty(xs)
        throw(ArgumentError("argument must not be empty"))
    end

    x_max = first(xs)
    for x in xs
        if isa(x, AbstractFloat) && isfinite(x)
            x_max = x
            break
        end
    end

    for x in xs
        if isa(x, AbstractFloat) && isfinite(x) && x > x_max
            x_max = x
        end
    end
    return x_max
end

function scientific_precision_heuristic(xs::AbstractArray{<:AbstractFloat})
    ys = [x == 0.0 ? 0.0 : round(10.0 ^ (z = log10(abs(Float64(x))); z - floor(z)); sigdigits=15)
          for x in xs if isfinite(x)]
    return plain_precision_heuristic(ys) + 1
end


function showoff(xs::AbstractArray{<:AbstractFloat}, style=:auto)
    x_min = concrete_minimum(xs)
    x_max = concrete_maximum(xs)
    x_min = Float64(x_min)
    x_max = Float64(x_max)

    if !isfinite(x_min) || !isfinite(x_max)
        return invoke(showoff,Tuple{AbstractArray,Symbol},xs,:none)
    end

    if style == :auto
        if x_max != x_min && abs(log10(x_max - x_min)) > 4
            style = :scientific
        else
            style = :plain
        end
    end

    if style == :plain
        precision = plain_precision_heuristic(xs)
        return String[format_fixed(x, precision) for x in xs]
    elseif style == :scientific
        precision = scientific_precision_heuristic(xs)
        return String[format_fixed_scientific(x, precision, false)
                      for x in xs]
    elseif style == :engineering
        precision = scientific_precision_heuristic(xs)
        return String[format_fixed_scientific(x, precision, true)
                      for x in xs]
    else
        throw(ArgumentError("$(style) is not a recongnized number format"))
    end
end

const superscript_numerals = ['⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹']


function showoff(ds::AbstractArray{T}, style=:none) where T<:Union{Date,DateTime}
    years = Set()
    months = Set()
    days = Set()
    hours = Set()
    minutes = Set()
    seconds = Set()
    for d in ds
        push!(years, Dates.year(d))
        push!(months, Dates.month(d))
        push!(days, Dates.day(d))
        push!(hours, Dates.hour(d))
        push!(minutes, Dates.minute(d))
        push!(seconds, Dates.second(d))
    end
    all_same_year         = length(years)   == 1
    all_one_month         = length(months)  == 1 && 1 in months
    all_one_day           = length(days)    == 1 && 1 in days
    all_zero_hour         = length(hours)   == 1 && 0 in hours
    all_zero_minute       = length(minutes) == 1 && 0 in minutes
    all_zero_seconds      = length(minutes) == 1 && 0 in minutes
    all_zero_milliseconds = length(minutes) == 1 && 0 in minutes

    # first label format
    label_months = false
    label_days = false
    f1 = "u d, yyyy"
    f2 = ""
    if !all_zero_seconds
        f2 = "HH:MM:SS.sss"
    elseif !all_zero_seconds
        f2 = "HH:MM:SS"
    elseif !all_zero_hour || !all_zero_minute
        f2 = "HH:MM"
    else
        if !all_one_day
            first_label_format = "u d yyyy"
        elseif !all_one_month
            first_label_format = "u yyyy"
        elseif !all_one_day
            first_label_format = "yyyy"
        end
    end
    if f2 != ""
        first_label_format = string(f1, " ", f2)
    else
        first_label_format = f1
    end

    labels = Vector{String}(undef, length(ds))
    labels[1] = Dates.format(ds[1], first_label_format)
    d_last = ds[1]
    for (i, d) in enumerate(ds[2:end])
        if Dates.year(d) != Dates.year(d_last)
            if all_one_day && all_one_month
                f1 = "yyyy"
            elseif all_one_day && !all_one_month
                f1 = "u yyyy"
            else
                f1 = "u d, yyyy"
            end
        elseif Dates.month(d) != Dates.month(d_last)
            f1 = all_one_day ? "u" : "u d"
        elseif Dates.day(d) != Dates.day(d_last)
            f1 = "d"
        else
            f1 = ""
        end

        if f2 != ""
            f = string(f1, " ", f2)
        elseif f1 != ""
            f = f1
        else
            f = first_label_format
        end

        labels[i+1] = Dates.format(d, f)
        d_last = d
    end

    return labels
end


end # module

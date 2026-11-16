import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
import ClimaUtilities.TimeManager: date
using Dates

"""
    AnnualCO2Series(csv_path; year_offset = 0, fixed_year = nothing)

Atmospheric CO2 (mol/mol) as a function of the simulation time, interpolated linearly
between annual means (ppm, columns `year,co2_ppm` of `csv_path`) placed at mid-year.

The nominal year is the forcing year plus `year_offset`, so that recycled forcing years
follow the CO2 of later nominal years. With `fixed_year`, the CO2 is that year's annual
mean at all times.
"""
struct AnnualCO2Series
    decimal_years::Vector{Float64}
    mol_per_mol::Vector{Float64}
    year_offset::Int
    fixed::Union{Nothing, Float64}
end

function AnnualCO2Series(csv_path; year_offset = 0, fixed_year = nothing)
    lines = filter(!isempty, readlines(csv_path))
    years = Int[]
    ppm = Float64[]
    for line in lines[2:end]
        fields = split(line, ',')
        push!(years, parse(Int, fields[1]))
        push!(ppm, parse(Float64, fields[2]))
    end
    fixed =
        isnothing(fixed_year) ? nothing :
        ppm[findfirst(==(fixed_year), years)] * 1e-6
    return AnnualCO2Series(years .+ 0.5, ppm .* 1e-6, year_offset, fixed)
end

decimal_year(d::DateTime) =
    Dates.year(d) +
    (Dates.value(d - DateTime(Dates.year(d))) / 1000) /
    (Dates.daysinyear(d) * 86400)

(series::AnnualCO2Series)(t) = co2_at(series, date(t))

"""
    co2_at(series::AnnualCO2Series, d::DateTime)

The CO2 (mol/mol) of `series` at the forcing date `d`.
"""
function co2_at(series::AnnualCO2Series, d::DateTime)
    isnothing(series.fixed) || return series.fixed
    x = decimal_year(d) + series.year_offset
    xs, ys = series.decimal_years, series.mol_per_mol
    x <= xs[1] && return ys[1]
    x >= xs[end] && return ys[end]
    i = searchsortedlast(xs, x)
    w = (x - xs[i]) / (xs[i + 1] - xs[i])
    return (1 - w) * ys[i] + w * ys[i + 1]
end

"""
    co2_time_varying_input(series::AnnualCO2Series)

The `TimeVaryingInput` of `series`, for the `c_co2` keyword of the forcing constructors.
"""
co2_time_varying_input(series::AnnualCO2Series) =
    TimeVaryingInput(t -> series(t))

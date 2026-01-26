module InverseFunctionsDatesExt

using Dates
import InverseFunctions: inverse

inverse(::typeof(Dates.datetime2epochms)) = Dates.epochms2datetime
inverse(::typeof(Dates.epochms2datetime)) = Dates.datetime2epochms
inverse(::typeof(Dates.date2epochdays)) = Dates.epochdays2date
inverse(::typeof(Dates.epochdays2date)) = Dates.date2epochdays

inverse(::typeof(datetime2unix)) = unix2datetime
inverse(::typeof(unix2datetime)) = datetime2unix
inverse(::typeof(datetime2julian)) = julian2datetime
inverse(::typeof(julian2datetime)) = datetime2julian

end

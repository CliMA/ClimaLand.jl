module StructUtilsTablesExt

using StructUtils, Tables

function StructUtils.applyeach(st::StructUtils.StructStyle, f, x::T) where {T <: Tables.AbstractRow}
    for nm in Tables.columnnames(x)
        ret = f(StructUtils.lowerkey(st, nm), StructUtils.lower(st, Tables.getcolumn(x, nm)))
        ret isa StructUtils.EarlyReturn && return ret
    end
    return StructUtils.defaultstate(st)
end

end # module
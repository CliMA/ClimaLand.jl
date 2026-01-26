# For a single and multi-file dataset
const NetCDFDataset =
    Union{NCDatasets.NCDataset, NCDatasets.CommonDataModel.MFDataset}

"""
    read_available_dates(ds::NCDatasets.NCDataset)

Return all the dates in the given NCDataset. The dates are read from the "time"
or "date" datasets. If none is available, return an empty vector.
"""
function read_available_dates(ds::NetCDFDataset)
    if "time" in keys(ds.dim)
        return Dates.DateTime.(
            reinterpret.(Ref(NCDatasets.DateTimeStandard), ds["time"][:]),
        )
    elseif "date" in keys(ds.dim)
        return Dates.DateTime.(string.(ds["date"][:]), Ref("yyyymmdd"))
    else
        return Dates.DateTime[]
    end
end

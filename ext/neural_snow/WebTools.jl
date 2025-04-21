using Downloads, DataFrames, Dates
export read_webpage_body, df_from_url

const extension_dir_path = @__DIR__()

"""
    read_webpage_body(url)

Utility function to grab the content of webpages and return in-memory as a String.
This function exists to circumnavigate weak dependencies on the HTTP package, though
this function subsequently has a more rigid scope of applicability. This returns a similar
ouput to `HTTP.get(url).body` as a String instead of including HTTP Response metadata or the
raw page content.

# Arguments
- `url::AbstractString`: the url for which to grab the page content.
"""
function read_webpage_body(url::AbstractString)
    @assert any(startswith.([url], ["http://", "https://"])) "Please provide a valid HTTP url.\n" #is there a better way to check the validity of the url without importing any more packages?
    #use the system temporary directory or write the temp file to here? aka, mktemp() or mktemp(extension_dir_path)?
    mktemp(extension_dir_path) do path, _
        try
            Downloads.request(url, output = path)
            if filesize(path) == 0
                error("Downloaded page at $(url) is empty or non-existent.")
            end
            body = String(read(path))  #like calling .body on a HTTP response
            return body
        catch e
            error("Download failed to fetch $(url)\n$(String(e))\n")
        end
    end
end

"""
    parse_cell(cell::AbstractString)

Utility function used in `df_from_url` to aid in the automated detection
of cell types in a web-based file being converted to a DataFrame. Supported
types currently include Ints, Float64s, Dates, and DateTimes, with the rest
being returned as Strings. Extend this function for additional automated
recognition capabilities.

# Arguments
- `cell::AbstractString`: the cell content from which to attempt automated type detection.
"""
function parse_cell(cell::AbstractString)
    isempty(cell) && return missing
    for T in (Int, Float64, Date, DateTime)
        val = tryparse(T, cell)
        val !== nothing && return val
    end
    return cell
end

"""
    df_from_url(url)

Utility function to grab the content of webpages and parse the page content
into a DataFrame. This is primarily used to handle the retrieval of data from
SNOTEL webpages within the DataTools module as part of the NeuralSnow extension to ClimaLand.
This function exists to circumnavigate weak dependencies on the HTTP and CSV packages,
executing similar functionality to `CSV.read(HTTP.get(url).body, DataFrame, delim = ",")`
with otherwise default arguments.

# Arguments
- `url::AbstractString`: the url for which to grab the page content.
- `comment::AbstractString`: the string for which commented-out lines are detected/ignored. Default is "#".
- `delim::AbstractString`: the string indicating the separation between entries. Default is ",".
"""
function df_from_url(url; comment = "#", delim = ",")
    body = read_webpage_body(url)
    lines = split(body, "\n", keepempty = false)
    data_rows = filter(
        l -> !(startswith(strip(l), comment) || isempty(strip(l))),
        lines,
    )
    header = split(strip(data_rows[1]), delim; keepempty = true)
    table =
        [split(strip(row), delim; keepempty = true) for row in data_rows[2:end]]
    columns = [getindex.(table, i) for i in 1:length(header)]
    parsed_cols = Array{Vector}(undef, length(header))
    for i in 1:length(header)
        parsed = parse_cell.(columns[i])
        T = eltype(skipmissing(parsed))
        if any(ismissing, parsed)
            parsed_cols[i] = convert(Vector{Union{Missing, T}}, parsed)
        else
            parsed_cols[i] = convert(Vector{T}, parsed)
        end
    end
    return DataFrame(parsed_cols, Symbol.(header))
end

"""Data utilities for running fluxnet site experiments. Provides the DataColumn
representation of a data column read from the FLUXNET data which includes 
a data status flag and units for the data. This enables the user to easily
convert units and check for missing data. Data is automatically filled when 
read in based on the quality control flags provided in the FLUXNET data or 
by simply filling with the mean value in the column."""

using Interpolations
using StatsBase

# Define the valid data statuses
@enum DataStatus complete = 1 absent = 2 incomplete = 3

"""
    replace_poor_quality_with_mean!(field, flag)

Replace values indicated to be poor quality by a fluxnet 
QC `flag` (Array) with the mean value
in the `field` (Array).

This uses the Fluxnet convention of 0 = measured, 1 = good quality
gap fill, 2 = medium quality, and 3 = poor quality, and replaces
data with QC flag = 2,3. 
"""
function replace_poor_quality_with_mean!(field, flag)
    good_indices = (flag .== 0) .|| (flag .== 1)
    fill_value = mean(field[good_indices])
    field[.~good_indices] .= fill_value
    return field
end

"""
    replace_replace_missing_with_mean_by_value!(field)

Replace values indicated to be missing
 with the mean value in the `field` (Array).

This uses the Fluxnet convention of a value of -9999 indicating
missing data.
"""
function replace_missing_with_mean_by_value!(field)
    good_indices = .~(field .== -9999)
    fill_value = mean(field[good_indices])
    field[.~good_indices] .= fill_value
    return field
end

"""
    column_status(driver_data::Matrix, column_name::String)

Checks that the `driver_data` matrix contains the column
specified by the string `column_name`; checks if that the column does not contain 
missing values. Returns the corresponding DataStatus (absent if column name 
is not present, incomplete if missing data is present, or complete if no missing
data is present).
"""
function column_status(driver_data::Matrix, column_name::String)
    col_dat = []
    column_names = driver_data[1, :]
    try
        col_dat = driver_data[2:end, column_names .== column_name]
    catch _
        return col_dat, absent
    end
    if isempty(col_dat)
        return col_dat, absent
    end
    try
        @assert(minimum(col_dat) > -9999)
    catch _
        if maximum(col_dat) == -9999
            return col_dat, absent
        end
        return col_dat, incomplete
    end
    return col_dat, complete
end

"""
    DataColumn

A struct for storing data columns along with thir units and 
a status flag which indicates the data is complete (but not neccessarily QC
checked), absent entirely, or incomplete.
"""
struct DataColumn
    "A nx1 matrix, where n is the number of data points, of data values"
    values::Matrix
    "The units of those values"
    units::String
    "The missing-data status of the column"
    status::DataStatus
end

"""
    filter_column(driver_data::Matrix, column_name::String, units::String)

Checks the `driver_data` for the present of `column_name`; returns a 
DataColumn struct with cleaned data and a complete status if present,
returns an empty DataColumn with status absent if not.

- Checks to see if a quality control flag 
is available for the data column, and if so, replaces poor quality values with
the mean value per the QC flag.
 - Checks for missing data and replaces with the mean value
- Returns an `absent` status if the column is not present.
"""
function filter_column(driver_data::Matrix, column_name::String, units::String)
    # Check that the data exists and read it in if so
    col_dat, status = column_status(driver_data, column_name)
    # Set missing data threshold above which column is treated as absent
    missing_threshold = 10.0
    # Set poor quality threshold above which the column undergoes no replacement using the quality flag
    quality_threshold = 10.0
    # if it does not exist, exit
    if status == absent
        @info "Warning: Data for $column_name is absent"
        return DataColumn(reshape([], 0, 0), "", status)
    elseif status == complete # all data is present but some may be poor quality
        # Check for a QC flag column and use it to clean data if possible
        column_QC = column_name * "_QC"
        QC_col, QCstatus = column_status(driver_data, column_QC)
        if QCstatus != absent
            # QC flag = [0,1] -> good data, [2,3] -> poorer data
            num_poor = count(x -> x >= 2, QC_col)
            percent_poor = 100.0 * num_poor / length(col_dat)
            # If the percent_poor is below a threshold, fill any missing data with the mean value per the QC flag, 
            # otherwise return data column with a specific warning but use poorer quality data directly
            if percent_poor < quality_threshold
                replace_poor_quality_with_mean!(col_dat, QC_col)
                @info "Warning: Data for $column_name $(round(percent_poor, sigdigits=3))% poor quality. Filled with mean value using QC flag"
                return DataColumn(col_dat, units, status)
            else
                @info "Warning: Data for $column_name $(round(percent_poor, sigdigits=3))% poor quality. Returning with no replacement."
                return DataColumn(col_dat, units, status)
            end
        else
            @info "Information: Data for $column_name is complete and no QC flag present"
            return DataColumn(col_dat, units, status)
        end
    elseif status == incomplete # data is missing
        # Replace values of -9999 with mean
        num_missing = count(x -> x == -9999, col_dat)
        percent_missing = 100.0 * num_missing / length(col_dat)
        if percent_missing > missing_threshold
            @info "Warning: Data for $column_name $(round(percent_missing,
                                        sigdigits=3))% has value of -9999. Treating as absent"
            return DataColumn(col_dat, "", absent)
        end
        replace_missing_with_mean_by_value!(col_dat)
        status = complete
        @info "Warning: Data for $column_name $(round(percent_missing,
                                sigdigits=3))% has value of -9999. Filled with mean value"
        return DataColumn(col_dat, units, status)
    end
end

"""
    transform_column(column::DataColumn,
                     transform::Function,
                     units::String,
                    )

Returns a new DataColumn struct with the values of the input `column`
transformed via the function `transform`. Useful for unit conversions. Sets 
the units of the column to the string units.
"""
function transform_column(
    column::DataColumn,
    transform::Function,
    units::String,
)
    return DataColumn(transform.(column.values), units, column.status)
end


"""
extract_variables(sv, variable_paths::Vector{String})

Extract multiple variables from saved simulation results.

Args:
    sv: Saved values from simulation 
    variable_paths: Vector of strings in format "module.variable" or "module.submodule.variable" 
                    (e.g., ["photosynthesis.GPP", "photosynthesis.OptVars.Vcmax25_opt"])

Returns:
    NamedTuple with variable names as keys and extracted time series as values
"""
function extract_variables(sv, variable_paths::Vector{String})
    results = NamedTuple()
    
    for path in variable_paths
        # Split the path into parts
        parts = split(path, ".")
        
        variable_data = []
        for k in 1:length(sv.saveval)
            try
                obj = sv.saveval[k]
                # Navigate through all parts of the path
                for part in parts
                    obj = getproperty(obj, Symbol(part))
                end
                push!(variable_data, parent(obj)[1])
            catch e
                # If extraction fails, push NaN
                push!(variable_data, NaN)
                println("Warning: Could not extract $path at timestep $k: $e")
            end
        end
        
        result_name = Symbol(parts[end])
        results = merge(results, NamedTuple{(result_name,)}((variable_data,)))
    end
    
    return results
end

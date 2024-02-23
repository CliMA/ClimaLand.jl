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


function replace_poor_by_mean!(field, flag, bounds::Vector)

    """
    this function is very similar to "replace_poor_quality_with_mean!" replaces field with flag values out of bounds range with average of field values whose flag are within the range
        # Arguments
        - `field`: An array of numerical data to be modified in place.
        - `flag`: An array of the same length as `field`, containing quality flags for each corresponding value in `field`.
        - `bounds`: A vector of two element where the first element is the lower bound and the second element is the upper bound of the flag values considered to indicate good quality.

        # Output
        - The function modifies `field` in place, replacing poor quality values with the mean of good quality values, and returns the modified `field`.
    """
    lower_bound, upper_bound = bounds
    good_indices = (flag .>= lower_bound) .&& (flag .<= upper_bound)
    fill_value = round(mean(field[good_indices]), digits = 2)
    field[.~good_indices] .= fill_value
    return field
end

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
    missing_threshold = 0.5
    # Set poor quality threshold above which the column undergoes no replacement using the quality flag
    quality_threshold = 50 # if more than 50% of the data in a site is poor quality
    # if it does not exist, exit
    bounds = [0, 1]#values between 0 and 1 are considered to be good quality

    if status == absent
        @info "Warning: Data for $column_name is absent"
        return DataColumn(reshape([], 0, 0), "", status)
    elseif status == complete # all data is present but some may be poor quality
        # Check for a QC flag column and use it to clean data if possible
        if column_name .== "GPP_DT_VUT_REF"
            column_QC = "NEE_VUT_REF_QC"
        elseif column_name .== "LE_CORR"
            column_QC = "LE_F_MDS_QC"
        elseif column_name .== "H_CORR"
            column_QC = "H_F_MDS_QC"
        else
            column_QC = column_name * "_QC"
        end
        QC_col, QCstatus = column_status(driver_data, column_QC)
        time = driver_data[2:end, 1]
        if QCstatus != absent
            # QC flag = [0,1] -> good data, [2,3] -> poorer data
            num_poor = count(x -> x >= bounds[2], QC_col)
            percent_poor = 100.0 * num_poor / length(col_dat)
            if percent_poor < quality_threshold
                # Set thr based on the presence of "SWC" or "GPP" in col_name
                thr =
                    contains(column_name, "SWC") ||
                    contains(column_name, "GPP") ? 0.0 : -900.0
                replace_with_longterm_mean!(col_dat, QC_col, time, bounds, thr)
                #replace_poor_quality_with_mean!(col_dat, QC_col,bounds)#just to check everything is alright
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
        if column_name .== "GPP_DT_VUT_REF"
            column_QC = "NEE_VUT_REF_QC"
        elseif column_name .== "LE_CORR"
            column_QC = "LE_F_MDS_QC"
        elseif column_name .== "H_CORR"
            column_QC = "H_F_MDS_QC"
        else
            column_QC = column_name * "_QC"
        end
        # Replace values of -9999 with mean
        num_missing = count(x -> x == -9999, col_dat)
        percent_missing = 100.0 * num_missing / length(col_dat)
        if percent_missing > missing_threshold
            @info "Warning: Data for $column_name $(round(percent_missing,
                                        sigdigits=3))% has value of -9999. Treating as absent"
            return DataColumn(col_dat, "", absent)
        end
        QC_col, QCstatus = column_status(driver_data, column_QC)
        thr =
            contains(column_name, "SWC") || contains(column_name, "GPP") ? 0.0 :
            -900.0
        replace_with_longterm_mean!(col_dat, QC_col, time, bounds, thr)
        #replace_missing_with_mean_by_value!(col_dat)
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

function replace_with_longterm_mean!(
    field,
    flag,
    time,
    bounds::Vector,
    thr::Float64,
)
    """
    this function replaces field data element that are below thr or have flag quality out of bounds range replaced with long term mean of same hour of data acquired on daily, weekly, biweekly, monthly, or yearly 
    # input arguments
    - `field`: the array of drivers
    - `flag`: an array of the same length as `field`, containing quality flags for each corresponding value in `field`.
    - `time `: the array same length as `field` that contains timestamps in the format of string "yyyymmddHHMM"
    - `bounds`: a vector of two element where the first element is the lower bound and the second element is the upper bound of the flag values considered to indicate good quality.
    - `thr`:  a single value that fields value below it get replaced regardless of flag value (for SWC and GPP it is often 0 and for other fields -900)
    # output argument
    - the modified `field` in place
    """
    function check_and_replace(criteria, field, flag, formatted_time)
        matching_indices = findall(criteria, formatted_time)
        filtered_flags = flag[matching_indices]
        filtered_field = field[matching_indices]
        if any(flg -> flg >= bounds[1] && flg <= bounds[2], filtered_flags)
            # Assuming replace_poor_quality_with_mean! modifies `field` in place for matching indices
            # and requires correct definition elsewhere
            replace_poor_by_mean!(filtered_field, filtered_flags, bounds)
            field[matching_indices] = filtered_field
            return true
        end
        return false
    end

    function update_flags_out_of_bounds!(field, flag, thr, out_of_bounds_value)
        """
        this function updates flag values manually
        checks if field values are nan, missing or below thr, sets the flag to be bad quality
        """

        for i in eachindex(field)
            if (ismissing(field[i]) || isnan(field[i]) || field[i] < thr)
                flag[i] = out_of_bounds_value
            end
        end
    end

    out_of_bounds_value = bounds[end] + 900
    formatted_time = DateTime.(string.(time), "yyyymmddHHMM")
    update_flags_out_of_bounds!(field, flag, thr, out_of_bounds_value)
    valid_values = filter(x -> !(x === missing || isnan(x) || x < thr), field)
    #first checks data belonging to the same hour of the year across all years

    if isempty(valid_values)
        error("This field doesn't contain any valid values to be gap-filled.")
    else
        indices = findall(
            x -> (
                ismissing(field[x]) ||
                isnan(field[x]) ||
                (flag[x] < bounds[1] || flag[x] > bounds[2])
            ),
            eachindex(field),
        )
        # filter indices
        for idx in indices
            dt = formatted_time[idx]
            time_criteria = [
                x ->
                    month(x) == month(dt) &&
                        day(x) == day(dt) &&
                        hour(x) == hour(dt),#same hour of the same day across all years
                x ->
                    month(x) == month(dt) &&
                        week(x) == week(dt) &&
                        hour(x) == hour(dt),#same hour of the same week of the year across all years
                x -> month(x) == month(dt) && hour(x) == hour(dt),#same hour of the same month across all years
                x -> month(x) == month(dt) && day(x) == day(dt),#average across the same day  of the same month
                x -> month(x) == month(dt) && week(x) == week(dt),#average across the same week of the same month
                x -> month(x) == month(dt),#same month of the year across all years
            ]
            for criterion in time_criteria
                if check_and_replace(criterion, field, flag, formatted_time)
                    return  # Exit after successful replacement
                end
            end
        end
    end

end

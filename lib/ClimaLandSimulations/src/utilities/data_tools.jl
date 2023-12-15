"""Data utilities for running fluxnet site experiments. Provides the DataColumn
representation of a data column read from the FLUXNET data which includes 
a data status flag and units for the data. This enables the user to easily
convert units and check for missing data. Data is automatically filled when 
read in based on the quality control flags provided in the FLUXNET data or 
by simply filling with the mean value in the column."""

# Define the valid data statuses
@enum DataStatus complete = 1 absent = 2 incomplete = 3

function replace_missing_with_mean!(field, flag)
    """Replace values indicated to be missing by a QC flag with the mean value
    in the column"""
    good_indices = (flag .== 0) .|| (flag .== 1)
    fill_value = mean(field[good_indices])
    field[.~good_indices] .= fill_value
    return field
end

function replace_missing_with_mean_by_value!(field)
    """Replace missing values indicated by -9999 in the column data with the 
    mean value in the column"""
    good_indices = .~(field .== -9999)
    fill_value = mean(field[good_indices])
    field[.~good_indices] .= fill_value
    return field
end

function check_column(file_matrix::Matrix, column_name::String)
    """Checks that the data file parsed into the data matrix contains the column
    specified by the string column, and that the column does not contain 
    missing values. Returns the corresponding DataStatus."""
    col_dat = []
    column_names = file_matrix[1, :]
    try
        col_dat = file_matrix[2:end, column_names .== column_name]
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

struct DataColumn
    """A struct for storing data columns along with thir units and 
    a status flag which indicates the data is complete, missing, or 
    incomplete."""
    values::Matrix
    units::String
    status::DataStatus
end

function DataColumn(file_matrix::Matrix, column_name::String, units::String)
    """Returns a DataColumn struct containing the data from the column
    specified by the string column_name, along with the units of the data and 
    a status flag which indicates the data is complete, missing, or
    incomplete."""
    col_dat, status = check_column(file_matrix, column_name)
    # Patch missing values with mean
    if status == incomplete
        num_missing = count(x -> x == -9999, col_dat)
        percent_missing = 100.0 * num_missing / length(col_dat)
        if percent_missing > 90.0
            return DataColumn(col_dat, units, absent)
        end
        @info "Warning: Data for $column_name $(round(percent_missing,
        sigdigits=3))% missing. Filling with mean value"
        replace_missing_with_mean_by_value!(col_dat)
        status = complete
    end
    return DataColumn(col_dat, units, status)
end

function VerifiedColumn(file_matrix::Matrix, column_name::String, units::String)
    """Returns a DataColumn struct with units and status for a column read from 
    the fluxnet data file. Checks to see if a quality control flag 
    is available for the data column, and if so, replaces missing values with
    the mean value per the QC flag."""
    # Check that the data exists and read it in
    col_dat, status = check_column(file_matrix, column_name)
    if status == absent
        return DataColumn(reshape([], 0, 0), "", status)
    end
    # Check for a QC flag column
    column_QC = column_name * "_QC"
    QC_col, QCstatus = check_column(file_matrix, column_QC)
    # If no QC flag, then fill and return the data column as is
    if QCstatus != absent && minimum(QC_col) != 2
        # Fill any missing data with the mean value per the QC flag
        num_missing = count(x -> x > 1, QC_col)
        if num_missing > 0
            perc_missing = 100.0 * num_missing / length(col_dat)
            @info "Warning: Data for $column_name $(
            round(perc_missing, sigdigits=3))% missing. Filling with mean value"
            replace_missing_with_mean!(col_dat, QC_col)
        end
        return DataColumn(col_dat, units, complete)
    end
    if status == incomplete
        num_missing = count(x -> x == -9999, col_dat)
        percent_missing = 100.0 * num_missing / length(col_dat)
        if percent_missing > 90.0
            return DataColumn(col_dat, "", absent)
        end
        @info "Warning: Data for $column_name $(round(percent_missing,
        sigdigits=3))% missing. Filling with mean value"
        status = complete
    end
    return DataColumn(col_dat, units, status)
end

function transform_column(
    column::DataColumn,
    transform::Function,
    units::String,
)
    """Returns a new DataColumn struct with the values of the input column
    transformed via the function transform. Useful for unit conversions. Sets 
    the units of the column to the string units."""
    return DataColumn(transform.(column.values), units, column.status)
end

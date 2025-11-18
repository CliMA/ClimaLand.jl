using Dates
using JLD2
import ClimaComms
using ClimaLand

include(joinpath(@__DIR__, "calculate_CWD.jl"))
include(joinpath(@__DIR__, "api.jl"))

"""
Generate CWD field file for model_interface.jl to use during calibration.
This avoids 17 ensemble members all trying to calculate CWD simultaneously.
"""
function generate_cwd_field_file(;
    start_year::Int = 2015,
    stop_year::Int = 2023
)
    @info "Generating CWD field file" start_year stop_year
    
    start_date = DateTime(start_year, 1, 1)
    stop_date = DateTime(stop_year, 12, 31)
    
    # Output filename (matches what model_interface.jl expects)
    if start_year == stop_year
        cwd_file = joinpath(@__DIR__, "cwd_era5_$(start_year).jld2")
    else
        cwd_file = joinpath(@__DIR__, "cwd_era5_$(start_year)_$(stop_year).jld2")
    end
    
    if isfile(cwd_file)
        @info "CWD field file already exists" cwd_file
        return cwd_file
    end
    
    # Manually construct ERA5 file paths for monthly structure
    base_path = "/net/sampo/data1/era5/forty_yrs_era5_land_forcing_data"
    era5_paths = String[]
    
    for year in start_year:stop_year
        for month in 1:12
            month_str = lpad(month, 2, '0')
            file_path = joinpath(base_path, "era_5_$year", "era5_forcing_data_$(year)_$(month_str).nc")
            if isfile(file_path)
                push!(era5_paths, file_path)
            else
                @warn "Missing file" file_path
            end
        end
    end
    
    @info "Found $(length(era5_paths)) ERA5 monthly file(s)"
    
    # Calculate CWD using vectorized function
    lons_cwd, lats_cwd, CWD_spatial = calculate_cwd_from_era5(
        era5_paths,
        start_date,
        stop_date,
    )
    
    # Save to file
    JLD2.jldsave(
        cwd_file;
        lons = lons_cwd,
        lats = lats_cwd,
        CWD_spatial,
        start_year,
        stop_year,
        created = now(),
        description = "CWD field for model_interface.jl (averaged $start_year-$stop_year)"
    )
    
    @info "Saved CWD field file" cwd_file size_MB=filesize(cwd_file)/1e6
    
    return cwd_file
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    generate_cwd_field_file(start_year = 2015, stop_year = 2023)
end
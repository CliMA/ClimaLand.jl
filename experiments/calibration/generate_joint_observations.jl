using JLD2
using Dates
using Statistics
using LinearAlgebra: diagm
import EnsembleKalmanProcesses as EKP

# DO NOT include api.jl or run_uspac_calibration.jl - they load the entire ClimaLand model
# which causes OOM during Julia initialization. Instead, hardcode config values.

"""
    generate_joint_observations(;
        lhf_filepath,
        smap_filepath,
        output_filepath,
        nelements,
        sample_date_ranges,
        lhf_weight,
        sm_weight
    )

Create joint observation vector by loading pre-generated ERA5 observations and combining with SMAP.
This avoids re-processing the huge ERA5 dataset.
"""
function generate_joint_observations(;
    lhf_filepath::String = "experiments/calibration/land_observation_vector.jld2",
    smap_filepath::String = "experiments/calibration/land_observation_vector_SMAP.jld2",
    output_filepath::String = "experiments/calibration/land_observation_vector_joint.jld2",
    nelements::Tuple{Int,Int} = (101, 15),
    sample_date_ranges::Vector = [("2016-12-1", "2019-9-1")],
    lhf_weight::Float64 = 0.5,
    sm_weight::Float64 = 0.5
)
    @info "Generating joint observations (loading pre-generated ERA5 + SMAP)" lhf_filepath smap_filepath
    
    # Load SMAP data selectively (only samples, not metadata - save memory)
    if !isfile(smap_filepath)
        error("SMAP observation file not found: $smap_filepath")
    end
    
    @info "Loading SMAP samples"
    sm_samples = JLD2.jldopen(smap_filepath, "r") do f
        f["observation_samples"]
    end
    
    @info "Loaded SMAP" n_periods=length(sm_samples)
    
    # Load ERA5 observations (already pre-processed, 6.4MB file)
    if !isfile(lhf_filepath)
        error("ERA5 observation file not found: $lhf_filepath. Run generate_observations.jl first.")
    end
    
    @info "Loading pre-generated ERA5 observations"
    n_periods = length(sample_date_ranges)
    
    # Verify both have same number of periods
    if length(sm_samples) != n_periods
        error("Period mismatch: ERA5 has $n_periods periods, SMAP has $(length(sm_samples)) periods")
    end
    
    # Build joint observations for each period
    joint_obs_vector = map(1:n_periods) do period_idx
        @info "Processing period $period_idx/$(n_periods)"
        
        # Load ERA5 observation for this period
        lhf_obs = JLD2.jldopen(lhf_filepath, "r") do file
            file["observation_vector"][period_idx]
        end
        
        # Get SMAP samples for this period (raw Float64 vector)
        sm_period_samples = sm_samples[period_idx]
        
        # Extract ERA5 samples
        lhf_samples = vcat(lhf_obs.samples...)
        
        @info "Period data" period=period_idx n_lhf=length(lhf_samples) n_sm=length(sm_period_samples)
        
        # Calculate scaling factors
        n_lhf = length(lhf_samples)
        n_sm = length(sm_period_samples)
        
        var_lhf = var(lhf_samples)
        var_sm = n_sm > 1 ? var(sm_period_samples) : 1.0
        
        # Scale by sqrt(n) to normalize total contribution, then by sqrt(var) to normalize units
        scale_lhf = sqrt(lhf_weight / n_lhf) / sqrt(var_lhf)
        scale_sm = sqrt(sm_weight / n_sm) / sqrt(var_sm)
        
        @info "Scaling" scale_lhf scale_sm ratio=scale_sm/scale_lhf
        
        # Scale observations
        lhf_scaled = lhf_samples .* Float32(scale_lhf)
        sm_scaled = sm_period_samples .* Float32(scale_sm)
        
        # Concatenate scaled observations: [LH, SMAP]
        combined_samples = vcat(lhf_scaled, sm_scaled)
        
        @info "Scaled and combined" period=period_idx n_total=length(combined_samples)
        
        # DON'T create EKP.Observation here - too memory intensive
        # Just return the scaled samples
        (samples = combined_samples, n_lhf = n_lhf, n_sm = n_sm)
    end
    
    @info "Created joint sample data" n_periods=length(joint_obs_vector)
    
    # Calculate statistics using first period for reporting
    first_samples = joint_obs_vector[1].samples
    n_sm_1 = joint_obs_vector[1].n_sm
    n_lhf_1 = joint_obs_vector[1].n_lhf
    n_total_first = length(first_samples)
    
    # Estimate statistics (approximate since data is scaled)
    println("\n===== Joint Observation Summary =====")
    println("LH flux observations (approx): $n_lhf_1")
    println("SMAP observations:             $n_sm_1")
    println("Total per period:              $n_total_first")
    println("Weighting: $(Int(lhf_weight*100))% LH / $(Int(sm_weight*100))% SMAP")
    
    # Save RAW SCALED SAMPLES (not EKP.Observation - too large)
    # The observation_map.jl will create EKP.Observation on-the-fly when needed
    JLD2.jldsave(
        output_filepath;
        joint_samples = [x.samples for x in joint_obs_vector],  # Just the scaled sample vectors
        n_lhf_per_period = [x.n_lhf for x in joint_obs_vector],
        n_sm_per_period = [x.n_sm for x in joint_obs_vector],
        nelements,
        sample_date_ranges,
        short_names = ["lhf", "sm_surface"],
        lhf_weight,
        sm_weight,
        created = now(),
        description = "Joint scaled samples (not EKP.Observation - create on load)"
    )
    
    @info "Saved joint samples" filepath=output_filepath n_periods=length(joint_obs_vector)
    
    return joint_obs_vector
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    # Hardcode config to avoid loading ClimaLand (causes OOM)
    sample_date_ranges = [("2016-12-1", "2019-9-1")]
    nelements = (101, 15)
    
    generate_joint_observations(
        lhf_filepath = "experiments/calibration/land_observation_vector.jld2",
        smap_filepath = "experiments/calibration/land_observation_vector_SMAP.jld2",
        output_filepath = "experiments/calibration/land_observation_vector_joint.jld2",
        nelements = nelements,
        sample_date_ranges = sample_date_ranges,
        lhf_weight = 0.5,
        sm_weight = 0.5
    )
end
#!/usr/bin/env julia
# =============================================================================
# Full-model OWUS calibration (CWD-only; static in time) against tower LE using EKP
# Parameters (6 total):
#   θ = [γ_f0, γ_fCWD,  γ_a0, γ_aCWD,  γ_b0, γ_bCWD]
# Forward model:
#   Build LandModel with OWUSConductanceCWDStatic(Γ, CWD), run, get diurnal LE (24)
# Observations:
#   Tower LE diurnal (24). Replace loader if you have a direct accessor.
# =============================================================================

using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Canopy
using ClimaLand.Simulations
import ClimaLand.FluxnetSimulations as FluxnetSimulations
import ClimaLand.Parameters as LP
import ClimaLand.LandSimVis as LandSimVis
import ClimaDiagnostics
import ClimaUtilities.TimeManager: date
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
#using CairoMakie
#CairoMakie.activate!()
using Statistics
using Logging
import Random
using Dates
#using ClimaAnalysis, GeoMakie, Printf, StatsBase

# ---------------- Helpers ----------------

# This function runs the model and computes diurnal average of latent heat flux
function G(θ)
    simulation = model(θ)
    lhf = get_lhf(simulation)
    observation =
        Float64.(
            get_diurnal_average(
                lhf,
                simulation.start_date,
                simulation.start_date + Day(20),
                stop_date,
            ),
        )
    return observation
end;

# Helper function: Extract latent heat flux from simulation diagnostics
function get_lhf(simulation)
    return ClimaLand.Diagnostics.diagnostic_as_vectors(
        simulation.diagnostics[1].output_writer,
        "lhf_30m_average",
    )
end;

# Helper function: Compute diurnal average of a variable
function get_diurnal_average(var, start_date, spinup_date, stop_date)
    (times, data) = var
    model_dates = if times isa Vector{DateTime}
        times
    else
        date.(times)
    end
    spinup_idx = findfirst(spinup_date .<= model_dates)
    stop_idx = findlast(model_dates .< stop_date)
    model_dates = model_dates[spinup_idx:stop_idx]
    data = data[spinup_idx:stop_idx]

    hour_of_day = Hour.(model_dates)
    mean_by_hour = [mean(data[hour_of_day .== Hour(i)]) for i in 0:23]
    return mean_by_hour
end;

# ---------------- Configuration and site setup ----------------

# Set random seed for reproducibility and floating point precision
rng_seed = 1234
rng = Random.MersenneTwister(rng_seed)
const FT = Float32

# Initialize land parameters and site configuration.
toml_dict = LP.create_toml_dict(FT)
site_ID = "US-MOz"
site_ID_val = FluxnetSimulations.replace_hyphen(site_ID);

# Get site-specific information: location coordinates, time offset, and sensor height.
(; time_offset, lat, long) =
    FluxnetSimulations.get_location(FT, Val(site_ID_val))
(; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, Val(site_ID_val));

# Get simulation start and stop dates in UTC; these must be included in the forcing data range Here we calibrate with only two months of data.
start_date = DateTime(2010, 5, 1)
stop_date = DateTime(2010, 7, 1)
Δt = 450.0; # seconds


# ---------------- Domain and forcing setup ----------------

# Create a column domain representing a 2-meter deep soil column with 10 vertical layers.
zmin = FT(-2)  # 2m depth
zmax = FT(0)   # surface
domain = Column(; zlim = (zmin, zmax), nelements = 10, longlat = (long, lat));

# Load prescribed atmospheric and radiative forcing from FLUXNET data
forcing = FluxnetSimulations.prescribed_forcing_fluxnet(
    site_ID,
    lat,
    long,
    time_offset,
    atmos_h,
    start_date,
    toml_dict,
    FT,
);

# Get Leaf Area Index (LAI) data from MODIS satellite observations.
LAI = ClimaLand.Canopy.prescribed_lai_modis(
    domain.space.surface,
    start_date,
    stop_date,
);

# ---------------- Model Setup ---------------


# Parameterization
# θ = [γ_f0, γ_fCWD, γ_a0, γ_aCWD, γ_b0, γ_bCWD]
function unpack_params(θ::AbstractVector)
    @assert length(θ) == 6
    γf0, γfC, γa0, γaC, γb0, γbC = θ
    return (γα0=γf0, γαC=γfC, γa0=γa0, γaC=γaC, γb0=γb0, γbC=γbC)
end

# Forward operator: build LandModel with OWUSConductanceCWDStatic

# Create an integrated land model that couples canopy, snow, soil, and soil CO2 components. 
# This comprehensive model allows us to simulate the full land surface system and its interactions.

function model(θ::AbstractVector)

    # unpack θ -> (α0, α1, a0, a1, b0, b1)
    α0, α1, a0, a1, b0, b1 = unpack_params(θ)

    # --- construct OWUS conductance from a parameters struct ---

    Γ_tuple = (FT(α0), FT(α1), FT(a0), FT(a1), FT(b0), FT(b1))  # NTuple{6,FT}
    
    owus_pars = Canopy.OWUSCWDStaticParameters{FT}(;
        Γ      = Γ_tuple,          # 3×2 matrix (rows: α,a,b; cols: intercept, log1p(CWD))
        cwd_mm = FT(CWD_MM),    # scalar site covariate (mm)
        # gsw_max = FT(Inf),    # optional, keep default unless you want a cap
    )
    

    #md # Set up ground conditions and define which components to simulate prognostically
    ground = ClimaLand.PrognosticGroundConditions{FT}()
    prognostic_land_components = (:canopy, :snow, :soil, :soilco2)

    #md # Prepare canopy domain and forcing
    canopy_domain = ClimaLand.Domains.obtain_surface_domain(domain)
    canopy_forcing = (; forcing.atmos, forcing.radiation, ground);

    #md # Set up photosynthesis parameters using the Farquhar model
    Vcmax25 = FT(Vcmax25)
    photosyn_defaults = Canopy.clm_photosynthesis_parameters(
        canopy_domain.space.surface,
    )
    photosynthesis = Canopy.FarquharModel{FT}(
        canopy_domain;
        photosynthesis_parameters = (;
            is_c3 = photosyn_defaults.is_c3,
            Vcmax25,
        ),
    )

    #conductance = Canopy.MedlynConductanceModel{FT}(canopy_domain; g1 = FT(141))
    conductance = Canopy.OWUSConductanceCWDStatic{FT}(owus_pars)

    #md # Create canopy model
    canopy = ClimaLand.Canopy.CanopyModel{FT}(
        canopy_domain,
        canopy_forcing,
        LAI,
        toml_dict;
        photosynthesis,
        prognostic_land_components,
        conductance,
    )

    #md # Create integrated land model
    land_model =
        LandModel{FT}(forcing, LAI, toml_dict, domain, Δt; canopy)

    #md # Set initial conditions from FLUXNET data
    set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
        site_ID,
        start_date,
        time_offset,
        land_model,
    )

    #md # Configure diagnostics to output latent heat fluxes half-hourly
    #md # Since fluxnet data is recorded every half hour, with the average of the half hour,
    #md # we need to be consistent here.
    output_vars = ["lhf"]
    diagnostics = ClimaLand.default_diagnostics(
        land_model,
        START_DATE;
        output_writer = ClimaDiagnostics.Writers.DictWriter(),
        output_vars,
        average_period = :halfhourly,
    )

    #md # Create and run the simulation
    simulation = Simulations.LandSimulation(
        start_date,
        stop_date,
        Δt,
        land_model;
        set_ic
        updateat = Second(Δt),
        user_callbacks = (),
        diagnostics,
    )
    solve!(simulation)
    return simulation
end

# ---------------- Experiment Setup ----------------

dataset = FluxnetSimulations.get_comparison_data(site_ID, time_offset)
observations = get_diurnal_average(
    (dataset.UTC_datetime, dataset.lhf),
    start_date,
    start_date + Day(20),
    stop_date,
);

noise_covariance = 0.05 * EKP.I;

# ---------------- Prior Distribution and Calibration Config ----------------

# Set up the prior distribution for the parameter and configure the ensemble Kalman inversion:
names = ["γ_f0","γ_fCWD","γ_a0","γ_aCWD","γ_b0","γ_bCWD"]
priors = PD.ParameterDistribution[
    PD.constrained_gaussian(names[i], 0.0, 0.30, -Inf, Inf) for i in 1:6
]
prior = PD.combine_distributions(priors)

# Set the ensemble size and number of iterations
ensemble_size = 10
N_iterations = 4;

# ---------------- Ensemble Kalman Inversion ----------------

initial_ensemble = EKP.construct_initial_ensemble(rng, prior, ensemble_size)
ensemble_kalman_process = EKP.EnsembleKalmanProcess(
    initial_ensemble, observations, noise_covariance, EKP.Inversion();
    scheduler = EKP.DataMisfitController(terminate_at = Inf, on_terminate = "continue"),
    rng,
)

Logging.with_logger(SimpleLogger(devnull, Logging.Error)) do
    for i in 1:N_iterations
        println("Iteration $i")
        θs = EKP.get_ϕ_final(prior, ensemble_kalman_process)
        G_ens = hcat([G(θs[:, j]...) for j in 1:ensemble_size]...)
        EKP.update_ensemble!(ensemble_kalman_process, G_ens)
    end
end;


# ---------------- Results Analysis and Visualization ----------------


θ̂ = EKP.get_ϕ_mean_final(prior, ekp)
Γ̂ = unpack_params(θ̂)
ŷ = model(θ̂)
rmse = sqrt(mean((ŷ .- observations).^2))

println("\n=== RESULT: OWUS (CWD-only, static) EKP vs tower LE ===")
println("Site=$(SITE)  Period=$(START_DATE)→$(STOP_DATE)")
@printf "Γ̂ (α0=%.3f, αC=%.3f,  a0=%.3f, aC=%.3f,  b0=%.3f, bC=%.3f)\n" Γ̂.γα0,Γ̂.γαC,Γ̂.γa0,Γ̂.γaC,Γ̂.γb0,Γ̂.γbC
@printf "RMSE(diurnal LE) = %.2f W m^-2\n" rmse


# plot
for i in 1:dim_size
    EKP.Visualize.plot_ϕ_over_iters(
        fig[1, i],
        ensemble_kalman_process,
        prior,
        i,
    )
end

EKP.Visualize.plot_error_over_iters(
    fig[1, dim_size + 1],
    ensemble_kalman_process,
)
#CairoMakie.save("constrained_params_and_error.png", fig)
fig







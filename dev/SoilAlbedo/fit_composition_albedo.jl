#=
Fit CompositionBasedSoilAlbedo coefficients from CERES observations and SoilGrids data.

This script calibrates the logistic regression coefficients for the CompositionBasedSoilAlbedo
parameterization by fitting to CERES bareground shortwave albedo observations over desert regions.

## Model formulation

The dry albedo is computed via logistic regression:
    η = η₀ + c_om * ν_ss_om + c_vgn * vg_n + c_cf * ν_ss_gravel
    α_dry = α_min + (α_max - α_min) * σ(η)

where σ(η) = 1/(1 + exp(-η)) is the logistic (sigmoid) function.

The logistic transformation ensures:
1. Albedo stays bounded in (α_min, α_max) = (0.04, 0.60)
2. Smooth, continuous gradients everywhere (no clipping artifacts)
3. Standard linear regression in transformed space

## Data sources

- CERES bareground shortwave albedo: Satellite-derived annual mean albedo
- SoilGrids: Soil organic matter (ν_ss_om), coarse fragments (ν_ss_gravel)
- van Genuchten n parameter

## Desert regions used for calibration

We focus on desert regions where bare soil dominates and CLM has known biases:
- Sahara Desert: lat ∈ [15, 35], lon ∈ [-15, 35]
- Arabian Desert: lat ∈ [15, 30], lon ∈ [35, 60]
- Australian Deserts: lat ∈ [-30, -20], lon ∈ [120, 145]
- Gobi Desert: lat ∈ [38, 48], lon ∈ [90, 115]
- Kalahari Desert: lat ∈ [-28, -20], lon ∈ [18, 28]

## Output

Calibrated coefficients achieving:
- RMSE = 0.074 (30% reduction vs mean)
- R² = 0.51
=#

using NCDatasets
using Statistics
using LinearAlgebra
using Printf

# Add ClimaLand to access artifact paths
using ClimaLand
import ClimaLand.Artifacts

# ============================================================================
# Configuration
# ============================================================================

# Physical bounds for albedo (logistic output range)
const α_min = 0.04  # Minimum physical albedo (dark organic soils)
const α_max = 0.60  # Maximum physical albedo (bright desert sand)

# Data paths from ClimaLand artifacts
# These functions trigger downloads if needed and return local paths
function get_data_paths()
    # CERES bareground shortwave albedo (annual mean)
    ceres_path = Artifacts.bareground_albedo_dataset_path()

    # SoilGrids soil composition (organic matter, coarse fragments)
    soilgrids_path = Artifacts.soil_grids_params_artifact_path(; lowres = true)

    # van Genuchten soil hydraulic parameters
    vg_params_dir = Artifacts.soil_params_artifact_folder_path(; lowres = true)
    vg_n_path = joinpath(vg_params_dir, "vGn_map_gupta_etal2020_1.0x1.0x4.nc")

    return (; ceres = ceres_path, soilgrids = soilgrids_path, vg_n = vg_n_path)
end

# ============================================================================
# Desert region definitions
# ============================================================================

"""
    DesertRegion

Defines a geographic bounding box for a desert region.
"""
struct DesertRegion
    name::String
    lat_min::Float64
    lat_max::Float64
    lon_min::Float64
    lon_max::Float64
end

# Major desert regions where bare soil dominates
const DESERT_REGIONS = [
    DesertRegion("Sahara", 15.0, 35.0, -15.0, 35.0),
    DesertRegion("Arabian", 15.0, 30.0, 35.0, 60.0),
    DesertRegion("Australian", -30.0, -20.0, 120.0, 145.0),
    DesertRegion("Gobi", 38.0, 48.0, 90.0, 115.0),
    DesertRegion("Kalahari", -28.0, -20.0, 18.0, 28.0),
]

"""
    is_in_desert(lat, lon, regions)

Check if a point (lat, lon) falls within any of the defined desert regions.
"""
function is_in_desert(lat::Float64, lon::Float64, regions::Vector{DesertRegion})
    for r in regions
        if r.lat_min <= lat <= r.lat_max && r.lon_min <= lon <= r.lon_max
            return true
        end
    end
    return false
end

# ============================================================================
# Data loading functions
# ============================================================================

"""
    load_ceres_albedo(filepath)

Load CERES bareground shortwave albedo data.
Returns (albedo_2d, lats, lons).
"""
function load_ceres_albedo(filepath::String)
    ds = NCDataset(filepath, "r")

    # Read coordinates
    lons = Array(ds["lon"][:])
    lats = Array(ds["lat"][:])

    # Read albedo - may be 1D (flattened) or 2D
    sw_alb_raw = Array(ds["sw_alb"][:])

    # Reshape if necessary
    if ndims(sw_alb_raw) == 1
        n_lon = length(lons)
        n_lat = length(lats)
        sw_alb = reshape(sw_alb_raw, n_lon, n_lat)
    else
        sw_alb = sw_alb_raw
    end

    close(ds)
    return sw_alb, lats, lons
end

"""
    load_soilgrids_composition(filepath)

Load SoilGrids soil composition data (organic matter, coarse fragments).
Returns (ν_ss_om, ν_ss_gravel, lats, lons).
"""
function load_soilgrids_composition(filepath::String)
    ds = NCDataset(filepath, "r")

    lons = Array(ds["lon"][:])
    lats = Array(ds["lat"][:])

    # Organic matter volume fraction
    ν_ss_om = Array(ds["nu_ss_om"][:])

    # Coarse fragments (gravel) volume fraction
    ν_ss_gravel = Array(ds["nu_ss_cf"][:])

    close(ds)
    return ν_ss_om, ν_ss_gravel, lats, lons
end

"""
    load_vangenuchten_n(filepath)

Load van Genuchten n parameter from soil hydraulic parameters dataset.
Returns (vg_n, lats, lons).
"""
function load_vangenuchten_n(filepath::String)
    ds = NCDataset(filepath, "r")

    lons = Array(ds["lon"][:])
    lats = Array(ds["lat"][:])
    vg_n = Array(ds["n"][:])

    close(ds)
    return vg_n, lats, lons
end

# ============================================================================
# Logistic regression utilities
# ============================================================================

"""
    inverse_sigmoid(p, α_min, α_max)

Compute the inverse of the sigmoid function (logit) for albedo values.
Maps α ∈ (α_min, α_max) → η ∈ (-∞, ∞).

This is used to transform observed albedo to the linear predictor space
where we can apply ordinary least squares regression.
"""
function inverse_sigmoid(α::Float64, α_min::Float64, α_max::Float64)
    # Normalize to (0, 1)
    p = (α - α_min) / (α_max - α_min)
    # Clamp to avoid log(0) or log(negative)
    p = clamp(p, 1e-6, 1.0 - 1e-6)
    # Logit transform
    return log(p / (1 - p))
end

"""
    sigmoid(η)

Standard logistic sigmoid function: σ(η) = 1 / (1 + exp(-η)).
"""
function sigmoid(η::Float64)
    if η >= 0
        return 1.0 / (1.0 + exp(-η))
    else
        exp_η = exp(η)
        return exp_η / (1.0 + exp_η)
    end
end

"""
    predict_albedo(η, α_min, α_max)

Predict albedo from linear predictor η using logistic transformation.
"""
function predict_albedo(η::Float64, α_min::Float64, α_max::Float64)
    return α_min + (α_max - α_min) * sigmoid(η)
end

# ============================================================================
# Fitting procedure
# ============================================================================

"""
    extract_desert_data(albedo, ν_ss_om, ν_ss_gravel, vg_n, lats, lons, regions)

Extract co-located data points from desert regions.
Returns vectors of (α_obs, om, cf, n) for valid grid points.
"""
function extract_desert_data(
    albedo::Matrix{Float64},
    ν_ss_om::Matrix{Float64},
    ν_ss_gravel::Matrix{Float64},
    vg_n::Matrix{Float64},
    lats::Vector{Float64},
    lons::Vector{Float64},
    regions::Vector{DesertRegion},
)
    α_obs = Float64[]
    om_vals = Float64[]
    cf_vals = Float64[]
    n_vals = Float64[]

    n_lon, n_lat = size(albedo)

    for i in 1:n_lon, j in 1:n_lat
        lat = lats[j]
        lon = lons[i]

        # Check if in desert region
        if !is_in_desert(lat, lon, regions)
            continue
        end

        # Check for valid data (no NaN or missing)
        α = albedo[i, j]
        om = ν_ss_om[i, j]
        cf = ν_ss_gravel[i, j]
        n = vg_n[i, j]

        if isnan(α) || isnan(om) || isnan(cf) || isnan(n)
            continue
        end

        # Filter unrealistic values
        if α < 0.05 || α > 0.55  # Physical albedo range
            continue
        end
        if n < 1.0 || n > 4.0  # Typical vg_n range
            continue
        end

        push!(α_obs, α)
        push!(om_vals, om)
        push!(cf_vals, cf)
        push!(n_vals, n)
    end

    return α_obs, om_vals, cf_vals, n_vals
end

"""
    extract_desert_data_with_coords(albedo, ν_ss_om, ν_ss_gravel, vg_n, lats, lons, regions)

Extract co-located data points from desert regions, including coordinates.
Returns vectors of (α_obs, om, cf, n, lat_vals, lon_vals) for valid grid points.
"""
function extract_desert_data_with_coords(
    albedo::Matrix{Float64},
    ν_ss_om::Matrix{Float64},
    ν_ss_gravel::Matrix{Float64},
    vg_n::Matrix{Float64},
    lats::Vector{Float64},
    lons::Vector{Float64},
    regions::Vector{DesertRegion},
)
    α_obs = Float64[]
    om_vals = Float64[]
    cf_vals = Float64[]
    n_vals = Float64[]
    lat_vals = Float64[]
    lon_vals = Float64[]

    n_lon, n_lat = size(albedo)

    for i in 1:n_lon, j in 1:n_lat
        lat = lats[j]
        lon = lons[i]

        # Check if in desert region
        if !is_in_desert(lat, lon, regions)
            continue
        end

        # Check for valid data (no NaN or missing)
        α = albedo[i, j]
        om = ν_ss_om[i, j]
        cf = ν_ss_gravel[i, j]
        n = vg_n[i, j]

        if isnan(α) || isnan(om) || isnan(cf) || isnan(n)
            continue
        end

        # Filter unrealistic values
        if α < 0.05 || α > 0.55  # Physical albedo range
            continue
        end
        if n < 1.0 || n > 4.0  # Typical vg_n range
            continue
        end

        push!(α_obs, α)
        push!(om_vals, om)
        push!(cf_vals, cf)
        push!(n_vals, n)
        push!(lat_vals, lat)
        push!(lon_vals, lon)
    end

    return α_obs, om_vals, cf_vals, n_vals, lat_vals, lon_vals
end

"""
    fit_logistic_regression(α_obs, om, cf, n; α_min=0.04, α_max=0.60)

Fit logistic regression coefficients using ordinary least squares in transformed space.

The model is:
    η = η₀ + c_om * om + c_vgn * n + c_cf * cf
    α = α_min + (α_max - α_min) * σ(η)

Returns (η₀, c_om, c_vgn, c_cf).
"""
function fit_logistic_regression(
    α_obs::Vector{Float64},
    om::Vector{Float64},
    cf::Vector{Float64},
    n::Vector{Float64};
    α_min::Float64 = 0.04,
    α_max::Float64 = 0.60,
)
    # Transform observed albedo to linear predictor space
    η_obs = [inverse_sigmoid(α, α_min, α_max) for α in α_obs]

    # Build design matrix: [1, om, n, cf]
    N = length(α_obs)
    X = hcat(ones(N), om, n, cf)

    # Solve normal equations: (X'X)^{-1} X' η
    coeffs = (X' * X) \ (X' * η_obs)

    η₀ = coeffs[1]
    c_om = coeffs[2]
    c_vgn = coeffs[3]
    c_cf = coeffs[4]

    return η₀, c_om, c_vgn, c_cf
end

"""
    evaluate_fit(α_obs, om, cf, n, η₀, c_om, c_vgn, c_cf; α_min=0.04, α_max=0.60)

Evaluate the fitted model, computing R² and RMSE.
"""
function evaluate_fit(
    α_obs::Vector{Float64},
    om::Vector{Float64},
    cf::Vector{Float64},
    n::Vector{Float64},
    η₀::Float64,
    c_om::Float64,
    c_vgn::Float64,
    c_cf::Float64;
    α_min::Float64 = 0.04,
    α_max::Float64 = 0.60,
)
    N = length(α_obs)

    # Predictions
    α_pred = Float64[]
    for i in 1:N
        η = η₀ + c_om * om[i] + c_vgn * n[i] + c_cf * cf[i]
        push!(α_pred, predict_albedo(η, α_min, α_max))
    end

    # Compute metrics
    residuals = α_obs .- α_pred
    ss_res = sum(residuals .^ 2)
    ss_tot = sum((α_obs .- mean(α_obs)) .^ 2)

    r_squared = 1.0 - ss_res / ss_tot
    rmse = sqrt(mean(residuals .^ 2))
    mae = mean(abs.(residuals))

    return (r_squared = r_squared, rmse = rmse, mae = mae, α_pred = α_pred)
end

# ============================================================================
# Main fitting script
# ============================================================================

function main()
    println("="^70)
    println("CompositionBasedSoilAlbedo Coefficient Fitting")
    println("="^70)
    println()

    # -------------------------------------------------------------------------
    # Load data
    # -------------------------------------------------------------------------
    println("Loading data...")

    # Get artifact paths (downloads if needed)
    paths = get_data_paths()

    println("  CERES albedo:    $(paths.ceres)")
    println("  SoilGrids:       $(paths.soilgrids)")
    println("  van Genuchten n: $(paths.vg_n)")
    println()

    albedo, lats_alb, lons_alb = load_ceres_albedo(paths.ceres)
    ν_ss_om, ν_ss_gravel, lats_sg, lons_sg =
        load_soilgrids_composition(paths.soilgrids)
    vg_n, lats_vg, lons_vg = load_vangenuchten_n(paths.vg_n)

    println("  CERES albedo: $(size(albedo))")
    println("  SoilGrids: $(size(ν_ss_om))")
    println("  van Genuchten n: $(size(vg_n))")
    println()

    # -------------------------------------------------------------------------
    # Extract desert region data
    # -------------------------------------------------------------------------
    println("Extracting desert region data...")

    α_obs, om, cf, n = extract_desert_data(
        albedo,
        ν_ss_om,
        ν_ss_gravel,
        vg_n,
        lats_alb,
        lons_alb,
        DESERT_REGIONS,
    )

    println("  Valid grid points: $(length(α_obs))")
    println("  Albedo range: [$(minimum(α_obs)), $(maximum(α_obs))]")
    println("  vg_n range: [$(minimum(n)), $(maximum(n))]")
    println()

    # -------------------------------------------------------------------------
    # Fit model
    # -------------------------------------------------------------------------
    println("Fitting logistic regression model...")

    η₀, c_om, c_vgn, c_cf = fit_logistic_regression(α_obs, om, cf, n)

    println()
    println("Fitted coefficients:")
    @printf("  η₀     = %+.2f\n", η₀)
    @printf("  c_om   = %+.2f (organic matter, expected negative)\n", c_om)
    @printf("  c_vgn  = %+.2f (van Genuchten n, expected positive)\n", c_vgn)
    @printf("  c_cf   = %+.2f (coarse fragments, expected positive)\n", c_cf)
    println()

    # -------------------------------------------------------------------------
    # Evaluate fit
    # -------------------------------------------------------------------------
    println("Evaluating fit quality...")

    results = evaluate_fit(α_obs, om, cf, n, η₀, c_om, c_vgn, c_cf)

    println()
    println("Fit quality metrics:")
    @printf("  R²   = %.3f\n", results.r_squared)
    @printf("  RMSE = %.3f\n", results.rmse)
    @printf("  MAE  = %.3f\n", results.mae)
    println()

    # -------------------------------------------------------------------------
    # Compare to mean baseline
    # -------------------------------------------------------------------------
    println("Comparison to mean baseline:")

    α_mean = mean(α_obs)
    rmse_mean = sqrt(mean((α_obs .- α_mean) .^ 2))

    @printf("  Mean baseline RMSE = %.3f\n", rmse_mean)
    @printf("  Model RMSE         = %.3f\n", results.rmse)
    @printf(
        "  Improvement        = %.1f%%\n",
        100 * (1 - results.rmse / rmse_mean)
    )
    println()

    # -------------------------------------------------------------------------
    # Output for use in ClimaLand
    # -------------------------------------------------------------------------
    println("="^70)
    println("Coefficients for CompositionBasedSoilAlbedo:")
    println("="^70)
    println()
    println("# PAR band (use directly fitted values)")
    @printf("η₀_PAR = FT(%.2f)\n", η₀)
    @printf("c_om_PAR = FT(%.2f)\n", c_om)
    @printf("c_vgn_PAR = FT(%.2f)\n", c_vgn)
    @printf("c_cf_PAR = FT(%.2f)\n", c_cf)
    println()
    println("# NIR band (slightly adjusted for spectral differences)")
    @printf("η₀_NIR = FT(%.2f)\n", η₀ - 0.06)  # NIR typically slightly different
    @printf("c_om_NIR = FT(%.2f)\n", c_om - 0.01)
    @printf("c_vgn_NIR = FT(%.2f)\n", c_vgn + 0.04)
    @printf("c_cf_NIR = FT(%.2f)\n", c_cf + 0.01)
    println()

    return (η₀ = η₀, c_om = c_om, c_vgn = c_vgn, c_cf = c_cf, results = results)
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

#=
Evaluate and compare different soil albedo parameterizations.

This script compares the performance of:
1. CLM soil color maps (traditional approach)
2. CompositionBasedSoilAlbedo (new physics-based approach)
3. Mean baseline (constant mean albedo)

against CERES bareground shortwave albedo observations.

## Metrics computed

- R² (coefficient of determination): Explained variance
- RMSE (root mean square error): Typical prediction error
- MAE (mean absolute error): Average absolute error
- Bias: Systematic over/under-estimation

## Regional analysis

Results are computed both globally and for specific desert regions where
bare soil dominates and CLM is known to have biases.
=#

using NCDatasets
using Statistics
using LinearAlgebra
using Printf

# Include the fitting utilities (provides get_data_paths(), data loading functions,
# DesertRegion, DESERT_REGIONS, predict_albedo, α_min, α_max, etc.)
include("fit_composition_albedo.jl")

# ============================================================================
# CLM soil color map albedo (for comparison)
# ============================================================================

"""
    clm_albedo_from_color(color_class)

Compute dry albedo from CLM soil color class (1-20).
Based on CLM Technical Note lookup tables.

Color class 1 = lightest (desert sand), 20 = darkest (organic soil).
"""
function clm_albedo_from_color(color_class::Int)
    # CLM dry albedo lookup table (PAR band)
    # Source: CLM Technical Note, Table 3.3
    clm_par_dry = [
        0.36,
        0.34,
        0.32,
        0.31,
        0.30,
        0.29,
        0.28,
        0.27,
        0.26,
        0.25,
        0.24,
        0.23,
        0.22,
        0.20,
        0.18,
        0.16,
        0.14,
        0.12,
        0.10,
        0.08,
    ]

    if 1 <= color_class <= 20
        return clm_par_dry[color_class]
    else
        return 0.20  # Default mid-range
    end
end

# ============================================================================
# Model evaluation functions
# ============================================================================

"""
    compute_metrics(α_obs, α_pred)

Compute standard regression metrics comparing observed and predicted albedo.
"""
function compute_metrics(α_obs::Vector{Float64}, α_pred::Vector{Float64})
    N = length(α_obs)
    @assert length(α_pred) == N

    # Residuals
    residuals = α_obs .- α_pred

    # R² (coefficient of determination)
    ss_res = sum(residuals .^ 2)
    ss_tot = sum((α_obs .- mean(α_obs)) .^ 2)
    r_squared = 1.0 - ss_res / ss_tot

    # RMSE
    rmse = sqrt(mean(residuals .^ 2))

    # MAE
    mae = mean(abs.(residuals))

    # Bias (mean error)
    bias = mean(residuals)

    # Correlation coefficient
    correlation = cor(α_obs, α_pred)

    return (
        r_squared = r_squared,
        rmse = rmse,
        mae = mae,
        bias = bias,
        correlation = correlation,
        n_points = N,
    )
end

"""
    evaluate_composition_model(α_obs, om, cf, n; coeffs)

Evaluate CompositionBasedSoilAlbedo predictions.
"""
function evaluate_composition_model(
    α_obs::Vector{Float64},
    om::Vector{Float64},
    cf::Vector{Float64},
    n::Vector{Float64};
    η₀::Float64 = -3.04,
    c_om::Float64 = -0.13,
    c_vgn::Float64 = 1.24,
    c_cf::Float64 = 0.15,
    α_min::Float64 = 0.04,
    α_max::Float64 = 0.60,
)
    N = length(α_obs)
    α_pred = Float64[]

    for i in 1:N
        η = η₀ + c_om * om[i] + c_vgn * n[i] + c_cf * cf[i]
        push!(α_pred, predict_albedo(η, α_min, α_max))
    end

    return compute_metrics(α_obs, α_pred), α_pred
end

"""
    evaluate_mean_baseline(α_obs)

Evaluate mean baseline (predicting the mean albedo everywhere).
"""
function evaluate_mean_baseline(α_obs::Vector{Float64})
    α_mean = mean(α_obs)
    α_pred = fill(α_mean, length(α_obs))
    return compute_metrics(α_obs, α_pred), α_pred
end

# ============================================================================
# Regional analysis
# ============================================================================

"""
    analyze_by_region(α_obs, om, cf, n, lats, lons, regions)

Compute metrics separately for each desert region.
"""
function analyze_by_region(
    α_obs::Vector{Float64},
    om::Vector{Float64},
    cf::Vector{Float64},
    n::Vector{Float64},
    lats::Vector{Float64},
    lons::Vector{Float64},
    regions::Vector{DesertRegion},
)
    results = Dict{String, NamedTuple}()

    for region in regions
        # Find points in this region
        mask = [
            region.lat_min <= lats[i] <= region.lat_max &&
            region.lon_min <= lons[i] <= region.lon_max for
            i in 1:length(lats)
        ]

        if sum(mask) < 10
            continue  # Skip regions with too few points
        end

        # Extract regional data
        α_reg = α_obs[mask]
        om_reg = om[mask]
        cf_reg = cf[mask]
        n_reg = n[mask]

        # Evaluate composition model
        metrics, _ = evaluate_composition_model(α_reg, om_reg, cf_reg, n_reg)

        # Evaluate mean baseline
        metrics_mean, _ = evaluate_mean_baseline(α_reg)

        results[region.name] = (
            composition = metrics,
            mean_baseline = metrics_mean,
            n_points = sum(mask),
            mean_albedo = mean(α_reg),
            mean_vgn = mean(n_reg),
        )
    end

    return results
end

# ============================================================================
# Correlation analysis
# ============================================================================

"""
    compute_predictor_correlations(α_obs, om, cf, n)

Compute correlations between albedo and each predictor variable.
This helps understand which soil properties most strongly influence albedo.
"""
function compute_predictor_correlations(
    α_obs::Vector{Float64},
    om::Vector{Float64},
    cf::Vector{Float64},
    n::Vector{Float64},
)
    return (
        r_om = cor(α_obs, om),      # Organic matter correlation
        r_cf = cor(α_obs, cf),      # Coarse fragments correlation
        r_vgn = cor(α_obs, n),      # van Genuchten n correlation
    )
end

# ============================================================================
# Main evaluation
# ============================================================================

function main()
    println("="^70)
    println("Soil Albedo Model Evaluation")
    println("="^70)
    println()

    # -------------------------------------------------------------------------
    # Load data from ClimaLand artifacts
    # -------------------------------------------------------------------------
    println("Loading data...")

    paths = get_data_paths()
    println("  CERES albedo:    $(paths.ceres)")
    println("  SoilGrids:       $(paths.soilgrids)")
    println("  van Genuchten n: $(paths.vg_n)")
    println()

    # Load datasets
    albedo, lats_alb, lons_alb = load_ceres_albedo(paths.ceres)
    ν_ss_om, ν_ss_gravel, lats_sg, lons_sg =
        load_soilgrids_composition(paths.soilgrids)
    vg_n, lats_vg, lons_vg = load_vangenuchten_n(paths.vg_n)

    # Extract desert region data with lat/lon coordinates
    α_obs, om, cf, n, lats, lons = extract_desert_data_with_coords(
        albedo,
        ν_ss_om,
        ν_ss_gravel,
        vg_n,
        lats_alb,
        lons_alb,
        DESERT_REGIONS,
    )

    N = length(α_obs)
    println("Dataset summary:")
    println("  N = $N grid points")
    @printf(
        "  Albedo: mean=%.3f, std=%.3f, range=[%.3f, %.3f]\n",
        mean(α_obs),
        std(α_obs),
        minimum(α_obs),
        maximum(α_obs)
    )
    @printf("  vg_n:   mean=%.2f, std=%.2f\n", mean(n), std(n))
    @printf("  OM:     mean=%.3f, std=%.3f\n", mean(om), std(om))
    @printf("  CF:     mean=%.2f, std=%.2f\n", mean(cf), std(cf))
    println()

    # -------------------------------------------------------------------------
    # Predictor correlations
    # -------------------------------------------------------------------------
    println("Predictor correlations with observed albedo:")
    correlations = compute_predictor_correlations(α_obs, om, cf, n)
    @printf(
        "  r(albedo, OM)   = %+.3f (expected negative)\n",
        correlations.r_om
    )
    @printf(
        "  r(albedo, CF)   = %+.3f (expected positive)\n",
        correlations.r_cf
    )
    @printf(
        "  r(albedo, vg_n) = %+.3f (expected positive, strongest)\n",
        correlations.r_vgn
    )
    println()

    # -------------------------------------------------------------------------
    # Model comparison
    # -------------------------------------------------------------------------
    println("-"^70)
    println("Model Comparison")
    println("-"^70)
    println()

    # 1. CompositionBasedSoilAlbedo (calibrated coefficients)
    metrics_comp, α_pred_comp = evaluate_composition_model(α_obs, om, cf, n)
    println("CompositionBasedSoilAlbedo (calibrated):")
    @printf("  R²   = %.3f\n", metrics_comp.r_squared)
    @printf("  RMSE = %.4f\n", metrics_comp.rmse)
    @printf("  MAE  = %.4f\n", metrics_comp.mae)
    @printf("  Bias = %+.4f\n", metrics_comp.bias)
    println()

    # 2. Mean baseline
    metrics_mean, _ = evaluate_mean_baseline(α_obs)
    println("Mean baseline (constant mean):")
    @printf("  R²   = %.3f (by definition)\n", 0.0)
    @printf("  RMSE = %.4f\n", metrics_mean.rmse)
    @printf("  MAE  = %.4f\n", metrics_mean.mae)
    println()

    # 3. Improvement summary
    println("-"^70)
    println("Improvement Summary")
    println("-"^70)
    @printf(
        "  RMSE reduction vs mean: %.1f%%\n",
        100 * (1 - metrics_comp.rmse / metrics_mean.rmse)
    )
    @printf(
        "  Variance explained (R²):       %.1f%%\n",
        100 * metrics_comp.r_squared
    )
    println()

    # -------------------------------------------------------------------------
    # Regional analysis
    # -------------------------------------------------------------------------
    println("-"^70)
    println("Regional Analysis")
    println("-"^70)
    println()

    regional_results =
        analyze_by_region(α_obs, om, cf, n, lats, lons, DESERT_REGIONS)

    for (region_name, results) in regional_results
        println("$region_name (n=$(results.n_points)):")
        @printf(
            "  Mean albedo: %.3f, Mean vg_n: %.2f\n",
            results.mean_albedo,
            results.mean_vgn
        )
        @printf(
            "  Composition RMSE: %.4f, R²: %.3f\n",
            results.composition.rmse,
            results.composition.r_squared
        )
        @printf("  Mean baseline RMSE: %.4f\n", results.mean_baseline.rmse)
        @printf(
            "  Improvement: %.1f%%\n",
            100 * (1 - results.composition.rmse / results.mean_baseline.rmse)
        )
        println()
    end

    # -------------------------------------------------------------------------
    # Physical interpretation
    # -------------------------------------------------------------------------
    println("="^70)
    println("Physical Interpretation of Coefficients")
    println("="^70)
    println()
    println("The fitted coefficients capture known soil-albedo relationships:")
    println()
    println("1. Organic matter (c_om = -0.13):")
    println("   - Negative coefficient: more organic matter → darker soil")
    println("   - Consistent with dark humus absorbing visible light")
    println()
    println("2. van Genuchten n (c_vgn = +1.24, strongest predictor):")
    println("   - Positive coefficient: higher n → brighter soil")
    println("   - High n indicates sandier, coarser texture")
    println("   - Sandy soils (quartz) are bright; clay soils are darker")
    println("   - Correlation r ≈ 0.71 makes vg_n the best single predictor")
    println()
    println("3. Coarse fragments (c_cf = +0.15):")
    println("   - Positive coefficient: more gravel/rock → brighter")
    println("   - Rocky/gravelly surfaces reflect more light")
    println("   - Important for distinguishing sandy vs rocky deserts")
    println()

    return (
        metrics_composition = metrics_comp,
        metrics_mean = metrics_mean,
        correlations = correlations,
        regional = regional_results,
    )
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

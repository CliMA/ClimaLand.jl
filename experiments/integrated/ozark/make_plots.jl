# Make all sort of plots with data and model output
# 1. Time series (e.g., C fluxes, h2o fluxes, energy fluxes, met drivers)
# 2. Seasonal pattern (i.e., monthly average and std)
# 3. Diurnal pattern (i.e., hourly average and std)
# 4. Response curves (e.g., NEE vs. PAR with VPD color and SWC brightness...)  
# 5. Energy conservation (i.e., Rn - G vs. L + H)
# 6. Water budget (i.e., cumulative ET vs. P)
# 7. Data quality plots (e.g., NEE vs. u* by T and SWC bins)
# 8. Fingerprint plots (showing both seasonality and diurnal pattern)
# 9. Wavelet coherence
# to do in another script: animations

# TO DO:  
# Script for plot utilities, Ed started this with plot_utils.jl

using ClimaLSM 
climalsm_dir = pkgdir(ClimaLSM)
include(joinpath(climalsm_dir, "experiments", "integrated", "ozark", "ozark.jl"))
savedir = joinpath(climalsm_dir, "experiments", "integrated", "ozark/") 
using CairoMakie # Draw vector graphics to SVG or PDF. High quality plots! 

# Get model output (using Ed's code)
outputs = save_model_outputs() # This is an array or arrays


# 1. Time series

# create an empty figure
fig = Figure(resolution = (1000, 1000)) # note: do not load Plots.jl in this branch 

# create empty axis, with a specific layout
ax_C = Axis(fig[1, 1]) # C fluxes
ax_W = Axis(fig[2, 1]) # h2o fluxes
ax_P = Axis(fig[3, 1]) # Precip & soil moisture
ax_T = Axis(fig[4, 1]) # air, canopy, and soil temperature

# add plots into axis ax_C
p_GPP_d
p_GPP_m 












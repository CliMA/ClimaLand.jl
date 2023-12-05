#=
include("integrated/ozark/ozark.jl")
include("integrated/ozark/input_df.jl")
include("integrated/ozark/model_output_dataframe.jl") # should rename df to outputs or climalsm
=# 

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

# TO DO:
# make white or dark background figures
# publication style and presentation style (bigger font etc.)

using ClimaLSM 
climalsm_dir = pkgdir(ClimaLSM)
savedir = joinpath(climalsm_dir, "experiments", "integrated", "ozark/") 
using CairoMakie # Draw vector graphics to SVG or PDF. High quality plots! 
using LaTeXStrings # To have latex labels

# drivers will be in drivers from Ed PR

# 1. Time series

# create an empty figure
fig = Figure(resolution = (2000, 2000)) # note: do not load Plots.jl in this branch 
fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

# create empty axis, with a specific layout
ax_C = Axis(fig[1, 1], ylabel = L"\text{GPP} \, (\mu\text{mol m}^{-2} \, \text{s}^{-1})") # C fluxes
ax_W = Axis(fig[2, 1]) # h2o fluxes
ax_P = Axis(fig[3, 1]) # Precip & soil moisture
ax_T = Axis(fig[4, 1]) # air, canopy, and soil temperature

# for time series, Makie should allow DateTime type soon (but not yet)
using PlotUtils: optimize_ticks
dateticks = optimize_ticks(df.DateTime[1], df.DateTime[end])[1][2:end-1] # first and last are weirdly placed

# add plots into axis ax_C
p_GPP_d = scatter!(ax_C, datetime2unix.(df.DateTime), df.GPP .* 1e6, color = :blue)
p_GPP_m = scatter!(ax_C, datetime2unix.(inputs.DateTime[index_t_start:index_t_end]), inputs.GPP[index_t_start:index_t_end] .*1e6, color = :black) 

ax_C.xticks[] = (datetime2unix.(dateticks) , Dates.format.(dateticks, "mm/dd"))

save("timeseries.png", fig)


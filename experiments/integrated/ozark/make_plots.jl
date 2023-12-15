#=
include("integrated/ozark/ozark.jl")
include("integrated/ozark/input_df.jl")
include("integrated/ozark/model_output_dataframe.jl") 

# TO DO:
# rename df to outputs or climalsm
# change unit from SI to what we want (e.g., umol m-2 s-1)
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
ax_W = Axis(fig[2, 1], ylabel = "ET (mm)") # h2o fluxes
ax_SWOUT = Axis(fig[3, 1], ylabel = L"\text{SW OUT} \, (\text{W} \, \text{m}^{-2})") # shortwave out 
# ax_T = Axis(fig[4, 1]) # air, canopy, and soil temperature

# for time series, Makie should allow DateTime type soon (but not yet)
# so the 2 lines of code below are a trick to be able to use DateTime - will be removed later
using PlotUtils: optimize_ticks
dateticks = optimize_ticks(df.DateTime[1], df.DateTime[end])[1][2:end-1] # first and last are weirdly placed

# add plots into axis ax_C
p_GPP_d = scatter!(ax_C, datetime2unix.(df.DateTime), df.GPP .* 1e6, color = :blue)
p_GPP_m = scatter!(ax_C, datetime2unix.(inputs.DateTime[index_t_start:index_t_end]), inputs.GPP[index_t_start:index_t_end] .*1e6, color = :black) 

# ax_W
p_ET_d = scatter!(ax_W, datetime2unix.(df.DateTime), (df.soil_evap .* 1e3 .* 24 .* 3600) .+ (df.transpiration .* 1e3 .* 24 .* 3600), color = :blue) # not sure about units
p_ET_m = scatter!(ax_W, datetime2unix.(inputs.DateTime[index_t_start:index_t_end]), inputs.LE[index_t_start:index_t_end] ./ (LSMP.LH_v0(earth_param_set) * 1000) .* (1e3 * 24 * 3600), color = :black) # not sure units

# ax_SW_OUT
p_SWOUT_d = scatter!(ax_SWOUT, datetime2unix.(df.DateTime), df.SW_out, color = :blue)
p_SWOUT_m = scatter!(ax_SWOUT, datetime2unix.(inputs.DateTime[index_t_start:index_t_end]), FT.(inputs.SW_OUT[index_t_start:index_t_end]), color = :black) 

# xticks
ax_C.xticks[] = (datetime2unix.(dateticks) , Dates.format.(dateticks, "mm/dd"))
ax_W.xticks[] = (datetime2unix.(dateticks) , Dates.format.(dateticks, "mm/dd"))
ax_SWOUT.xticks[] = (datetime2unix.(dateticks) , Dates.format.(dateticks, "mm/dd"))

axislegend(ax_C, [p_GPP_d, p_GPP_m], ["Observations", "ClimaLSM"], "", position = :rt, orientation = :horizontal)

save("timeseries.svg", fig)





















# Fingerplot
fig = Figure(resolution = (2000, 2000))
fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

# create empty axis, with a specific layout
ax_C = Axis(fig[1, 1], ylabel = "Hour of the day", xlabel = "Date", title = L"\text{GPP} \, (\mu\text{mol m}^{-2} \, \text{s}^{-1})") # C fluxes
ax_W = Axis(fig[2, 1]) # h2o fluxes
ax_P = Axis(fig[3, 1]) # Precip & soil moisture
ax_T = Axis(fig[4, 1]) # air, canopy, and soil temperature

# for time series, Makie should allow DateTime type soon (but not yet)
# so the 2 lines of code below are a trick to be able to use DateTime - will be removed later
using PlotUtils: optimize_ticks
dateticks = optimize_ticks(inputs.DateTime[1], inputs.DateTime[end])[1][2:end-1] # first and last are weirdly placed

# Fingerprint plot
hm_GPP = heatmap!(ax_C, datetime2unix.(DateTime.(Date.(inputs.DateTime))), hour.(inputs.DateTime) .+ (minute.(inputs.DateTime) ./ 60), inputs.GPP .* 1e6)
Colorbar(fig[1, 2], hm_GPP)

ax_C.xticks[] = (datetime2unix.(dateticks) , Dates.format.(dateticks, "mm/dd"))

save("Fingerprint.svg", fig)






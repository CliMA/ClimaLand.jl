#= run these lines for test / development

ARGS = ["US-MOz"]
include("integrated/fluxnet/setup.jl")
include("integrated/fluxnet/run_fluxnet.jl")
include("integrated/fluxnet/inputs_dataframe.jl")
include("integrated/fluxnet/climalsm_output_dataframe.jl") 

# TO DO:
# change unit from SI to what we want (e.g., umol m-2 s-1)
# soil moisture is currently deepest I think, change the depth

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
savedir =
    joinpath(climalsm_dir, "experiments", "integrated", "fluxnet/figures/")
using CairoMakie # Draw vector graphics to SVG or PDF. High quality plots! 
using LaTeXStrings # To have latex labels
using PlotUtils: optimize_ticks

# drivers will be in drivers from Ed PR

# 1. Time series of GPP, ET and SW_OUT
function timeseries_fluxes_fig(inputs, climalsm) # will run for any inputs or climalsm output of FLUXNET sites
    # create an empty figure
    fig = Figure(size = (1000, 1000)) # note: do not load Plots.jl in this branch (it is loading in plot_utils)
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    # create empty axis, with a specific layout
    ax_C = Axis(
        fig[1, 1],
        ylabel = L"\text{GPP} \, (\mu\text{mol m}^{-2} \, \text{s}^{-1})",
    ) # C fluxes
    ax_W = Axis(fig[2, 1], ylabel = L"\text{ET (mm)}") # h2o fluxes
    ax_SWOUT = Axis(
        fig[3, 1],
        ylabel = L"\text{SW OUT} \, (\text{W} \, \text{m}^{-2})",
    ) # shortwave out 
    # ax_T = Axis(fig[4, 1]) # air, canopy, and soil temperature

    # for time series, Makie should allow DateTime type soon (but not yet)
    # so the 2 lines of code below are a trick to be able to use DateTime - will be removed later
    dateticks =
        optimize_ticks(climalsm.DateTime[1], climalsm.DateTime[end])[1][2:(end - 1)] # first and last are weirdly placed

    # add plots into axis ax_C
    p_GPP_m = CairoMakie.scatter!(
        ax_C,
        datetime2unix.(climalsm.DateTime),
        climalsm.GPP .* 1e6,
        color = :blue,
    )
    p_GPP_d = CairoMakie.scatter!(
        ax_C,
        datetime2unix.(inputs.DateTime[index_t_start:index_t_end]),
        inputs.GPP[index_t_start:index_t_end] .* 1e6,
        color = :black,
    )

    # ax_W
    p_ET_m = CairoMakie.scatter!(
        ax_W,
        datetime2unix.(climalsm.DateTime),
        (climalsm.vapor_flux .* 1e3 .* 24 .* 3600) .+
        (climalsm.transpiration .* 1e3 .* 24 .* 3600),
        color = :blue,
    ) # not sure about units
    p_ET_d = CairoMakie.scatter!(
        ax_W,
        datetime2unix.(inputs.DateTime[index_t_start:index_t_end]),
        inputs.LE[index_t_start:index_t_end] ./
        (LSMP.LH_v0(earth_param_set) * 1000) .* (1e3 * 24 * 3600),
        color = :black,
    ) # not sure units

    # ax_SW_OUT
    p_SWOUT_m = CairoMakie.scatter!(
        ax_SWOUT,
        datetime2unix.(climalsm.DateTime),
        climalsm.SW_out,
        color = :blue,
    )
    p_SWOUT_d = CairoMakie.scatter!(
        ax_SWOUT,
        datetime2unix.(inputs.DateTime[index_t_start:index_t_end]),
        FT.(inputs.SW_OUT[index_t_start:index_t_end]),
        color = :black,
    )

    # xticks
    ax_C.xticks[] =
        (datetime2unix.(dateticks), Dates.format.(dateticks, "mm/dd"))
    ax_W.xticks[] =
        (datetime2unix.(dateticks), Dates.format.(dateticks, "mm/dd"))
    ax_SWOUT.xticks[] =
        (datetime2unix.(dateticks), Dates.format.(dateticks, "mm/dd"))

    axislegend(
        ax_C,
        [p_GPP_d, p_GPP_m],
        ["Observations", "ClimaLSM"],
        "",
        position = :rt,
        orientation = :horizontal,
    )

    CairoMakie.ylims!(ax_C, (0, 40))
    CairoMakie.ylims!(ax_W, (0, 30))
    CairoMakie.ylims!(ax_SWOUT, (0, 200))

    [
        CairoMakie.xlims!(
            axes,
            (
                datetime2unix(climalsm.DateTime[1]),
                datetime2unix(climalsm.DateTime[end]),
            ),
        ) for axes in [ax_C, ax_W, ax_SWOUT]
    ]

    fig
    return fig
end

#= to test. These will be in another script though.
fig = timeseries_fluxes_fig(inputs, climalsm) 
save(joinpath(savedir, "timeseries_fluxes.pdf"), fig)
=#

# 2. Time series of SWC, Precip, moisture stress, stomatal conductance
function timeseries_H2O_fig(inputs, climalsm) # will run for any inputs or climalsm output of FLUXNET sites
    # create an empty figure
    fig = Figure(size = (1000, 1000)) # note: do not load Plots.jl in this branch (it is loading in plot_utils)
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    # create empty axis, with a specific layout
    ax_H2O =
        Axis(fig[1, 1], ylabel = L"\theta \, (\text{m}^{3} \, \text{m}^{-3})") # soil moisture
    ax_H2O_rain = Axis(
        fig[1, 1],
        ylabel = L"\text{Rainfall (mm)}",
        yaxisposition = :right,
        ygridvisible = false,
        ylabelcolor = :blue,
        yticklabelcolor = :blue,
    )
    ax_MS = Axis(fig[2, 1], ylabel = L"\text{Moisture stress}") # moisture stress
    ax_SC = Axis(
        fig[3, 1],
        ylabel = L"\text{Stomatal conductance} \, (\text{mol m}^{-2} \, \text{s}^{-1})",
    ) # stomatal conductance 
    # ax_T = Axis(fig[4, 1]) # air, canopy, and soil temperature

    # for time series, Makie should allow DateTime type soon (but not yet)
    # so the 2 lines of code below are a trick to be able to use DateTime - will be removed later
    dateticks =
        optimize_ticks(climalsm.DateTime[1], climalsm.DateTime[end])[1][2:(end - 1)] # first and last are weirdly placed

    # add plots into axis ax_H2O
    p_H2O_m = CairoMakie.scatter!(
        ax_H2O,
        datetime2unix.(climalsm.DateTime),
        climalsm.θ_l,
        color = :green,
    )
    p_H2O_d = CairoMakie.scatter!(
        ax_H2O,
        datetime2unix.(inputs.DateTime[index_t_start:index_t_end]),
        inputs.SWC[index_t_start:index_t_end],
        color = :black,
    )

    # Rain on secondary axis
    p_rain = barplot!(
        ax_H2O_rain,
        datetime2unix.(inputs.DateTime[index_t_start:index_t_end]),
        inputs.P[index_t_start:index_t_end],
        color = :blue,
    )

    # Moisture stress
    p_MS = CairoMakie.scatter!(
        ax_MS,
        datetime2unix.(climalsm.DateTime),
        climalsm.β,
        color = :green,
    ) # not sure about units

    # Stomatal conductance
    p_SC = CairoMakie.scatter!(
        ax_SC,
        datetime2unix.(climalsm.DateTime),
        climalsm.gs,
        color = :green,
    )

    # xticks
    ax_H2O.xticks[] =
        (datetime2unix.(dateticks), Dates.format.(dateticks, "mm/dd"))
    ax_H2O_rain.xticks[] =
        (datetime2unix.(dateticks), Dates.format.(dateticks, "mm/dd"))
    hidexdecorations!(ax_H2O_rain)
    ax_MS.xticks[] =
        (datetime2unix.(dateticks), Dates.format.(dateticks, "mm/dd"))
    ax_SC.xticks[] =
        (datetime2unix.(dateticks), Dates.format.(dateticks, "mm/dd"))

    axislegend(
        ax_H2O,
        [p_H2O_d, p_H2O_m],
        ["Observations", "ClimaLSM"],
        "",
        position = :rt,
        orientation = :horizontal,
    )

    #CairoMakie.ylims!(ax_C, (0, 40))
    #CairoMakie.ylims!(ax_W, (0, 30))
    #CairoMakie.ylims!(ax_SWOUT, (0, 200))

    [
        CairoMakie.xlims!(
            axes,
            (
                datetime2unix(climalsm.DateTime[1]),
                datetime2unix(climalsm.DateTime[end]),
            ),
        ) for axes in [ax_H2O, ax_H2O_rain, ax_MS, ax_SC]
    ]

    fig
    return fig
end

#= to test. These will be in another script though.
fig = timeseries_H2O_fig(inputs, climalsm) 
save(joinpath(savedir, "timeseries_H2O.pdf"), fig)
=#

# 3. Fingerprint plot
function fingerprint_fig(inputs, climalsm)
    fig = Figure(size = (1000, 1000))
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    # create empty axis, with a specific layout
    ax_C = Axis(
        fig[1, 1],
        ylabel = "Hour of the day",
        xlabel = "Date",
        title = L"\text{GPP} \, (\mu\text{mol m}^{-2} \, \text{s}^{-1})",
    ) # C fluxes
    ax_W = Axis(fig[2, 1]) # h2o fluxes
    ax_P = Axis(fig[3, 1]) # Precip & soil moisture
    ax_T = Axis(fig[4, 1]) # air, canopy, and soil temperature

    # for time series, Makie should allow DateTime type soon (but not yet)
    # so the 2 lines of code below are a trick to be able to use DateTime - will be removed later
    dateticks =
        optimize_ticks(inputs.DateTime[1], inputs.DateTime[end])[1][2:(end - 1)] # first and last are weirdly placed

    # Fingerprint plot
    hm_GPP = CairoMakie.heatmap!(
        ax_C,
        datetime2unix.(DateTime.(Date.(inputs.DateTime))),
        hour.(inputs.DateTime) .+ (minute.(inputs.DateTime) ./ 60),
        inputs.GPP .* 1e6,
    )
    Colorbar(fig[1, 2], hm_GPP)

    ax_C.xticks[] =
        (datetime2unix.(dateticks), Dates.format.(dateticks, "mm/dd"))

    fig
    return fig
end

#=
fig = fingerprint_fig(inputs, climalsm)
save(joinpath(savedir, "fingerprint.pdf"), fig)
=#

# 4. Diurnals, with quantiles, for C, h2o, energy

# Energy: H, L, G, LW_OUT, SW_OUT
# inputs.G, inputs.H, inputs.LE, inputs.SW_OUT 
# climalsm.lhf, climalsm.shf, climalsm.LW_out, climalsm.SW_OUT # need G

# C: GPP, ER (AR + HR)
# inputs.RECO, inputs.GPP
# climalsm.Ra, climalsm.GPP # need Rh

# drivers: Temperature (soil, air, canopy)
# inputs.TS, inputs.TA
# climalsm.T # need canopy and soil T?



using Statistics

function diurnal(datetime, data)
    hourlyquantile = [
        quantile(data[hour.(datetime) .== h], [0.25, 0.5, 0.75]) for h in 0:1:23
    ]
    return hourlyquantile
end

function diurnal_plot!(
    fig,
    ax,
    datetime,
    data,
    color;
    alpha = 0.3,
    linestyle = :solid,
)
    fig
    hourlyquantile = diurnal(datetime, data)
    diurnal_p = CairoMakie.lines!(
        ax,
        0.5:1:23.5,
        getindex.(hourlyquantile[1:24], 2),
        color = color,
        linestyle = linestyle,
    )
    diurnal_q = CairoMakie.band!(
        ax,
        0.5:1:23.5,
        getindex.(hourlyquantile[1:24], 1),
        getindex.(hourlyquantile[1:24], 3),
        color = (color, alpha),
    )
    return diurnal_p
end

function diurnals_fig(inputs, climalsm)
    fig = Figure(size = (1000, 1000))
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    ax_C = Axis(
        fig[1, 1],
        ylabel = L"\text{CO}_{2} \, (\mu\text{mol m}^{-2} \, \text{s}^{-1})",
    ) # C fluxes
    ax_W = Axis(fig[2, 1], ylabel = L"\text{H}_{2}\text{O} \, \text{(mm)}") # h2o fluxes
    ax_E = Axis(
        fig[3, 1],
        ylabel = L"\text{Radiation} \, (\text{W} \, \text{m}^{-2})",
        xlabel = L"\text{Hour of the day}",
        xgridvisible = false,
    ) # shortwave out 

    # CO2 fluxes
    # model
    p_GPP_m =
        diurnal_plot!(fig, ax_C, climalsm.DateTime, climalsm.GPP .* 1e6, :green)
    diurnal_plot!(fig, ax_C, climalsm.DateTime, climalsm.Ra .* 1e6, :black)
    # data
    p_GPP_d = diurnal_plot!(
        fig,
        ax_C,
        inputs.DateTime[index_t_start:index_t_end],
        inputs.GPP[index_t_start:index_t_end] .* 1e6,
        :green,
        alpha = 0.1,
        linestyle = :dot,
    )

    # H2O fluxes
    # model
    p_ET_m = diurnal_plot!(
        fig,
        ax_W,
        climalsm.DateTime,
        climalsm.transpiration .* 1e3 .* 24 .* 3600,
        :blue,
    )
    # data
    p_ET_d = diurnal_plot!(
        fig,
        ax_W,
        inputs.DateTime[index_t_start:index_t_end],
        inputs.LE[index_t_start:index_t_end] ./
        (LSMP.LH_v0(earth_param_set) * 1000) .* (1e3 * 24 * 3600),
        :blue,
        alpha = 0.1,
        linestyle = :dot,
    )

    # Energy fluxes
    # model
    # diurnal_plot!(fig, ax_E, climalsm.DateTime, climalsm.LW_out, :red)
    p_SWout_m =
        diurnal_plot!(fig, ax_E, climalsm.DateTime, climalsm.SW_out, :red)
    # data
    p_SWout_d = diurnal_plot!(
        fig,
        ax_E,
        inputs.DateTime[index_t_start:index_t_end],
        FT.(inputs.SW_OUT[index_t_start:index_t_end]),
        :red,
        alpha = 0.1,
        linestyle = :dot,
    )

    [CairoMakie.xlims!(axes, (0, 24)) for axes in [ax_C, ax_W, ax_E]]

    axislegend(
        ax_C,
        [p_GPP_d, p_GPP_m],
        ["Observations", "ClimaLSM"],
        "",
        position = :rt,
        orientation = :horizontal,
    )
    axislegend(
        ax_W,
        [p_ET_d, p_ET_m],
        ["Observations", "ClimaLSM"],
        "",
        position = :rt,
        orientation = :horizontal,
    )
    axislegend(
        ax_E,
        [p_SWout_d, p_SWout_m],
        ["Observations", "ClimaLSM"],
        "",
        position = :rt,
        orientation = :horizontal,
    )

    hidexdecorations!(ax_C)
    hidexdecorations!(ax_W)

    fig
    return fig
end

#=
fig = diurnals_fig(inputs, climalsm)
save(joinpath(savedir, "diurnals.pdf"), fig)
=#

# 5. Cumulative P and ET

# 6. Energy balance closure (L + H = Rn - G)

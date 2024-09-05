export timeseries_fluxes_fig,
    timeseries_H2O_fig,
    fingerprint_fig,
    diurnal,
    diurnal_plot!,
    diurnals_fig,
    make_plots

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

# 1. Time series of GPP, ET and SW_OUT
function timeseries_fluxes_fig(
    inputs,
    climaland,
    earth_param_set;
    dashboard = false,
) # will run for any inputs or climaland output of FLUXNET sites

    if dashboard == true
        WGLMakie.activate!() # for dashboards
    end

    index_t_start = findfirst(isequal(climaland.DateTime[1]), inputs.DateTime)
    index_t_end = findfirst(isequal(climaland.DateTime[end]), inputs.DateTime)

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

    # for time series, CairoMakie should allow DateTime type soon (but not yet)
    # so the 2 lines of code below are a trick to be able to use DateTime - will be removed later
    dateticks =
        optimize_ticks(climaland.DateTime[1], climaland.DateTime[end])[1][2:(end - 1)] # first and last are weirdly placed

    # add plots into axis ax_C
    p_GPP_m = lines!(
        ax_C,
        datetime2unix.(climaland.DateTime),
        climaland.GPP .* 1e6,
        color = :blue,
    )
    p_GPP_d = lines!(
        ax_C,
        datetime2unix.(inputs.DateTime[index_t_start:index_t_end]),
        inputs.GPP[index_t_start:index_t_end] .* 1e6,
        color = :black,
    )

    # ax_W
    p_ET_m = lines!(
        ax_W,
        datetime2unix.(climaland.DateTime),
        (climaland.vapor_flux_liq .* 1e3 .* 24 .* 3600) .+
        (climaland.transpiration .* 1e3 .* 24 .* 3600),
        color = :blue,
    ) # not sure about units
    p_ET_d = lines!(
        ax_W,
        datetime2unix.(inputs.DateTime[index_t_start:index_t_end]),
        inputs.LE[index_t_start:index_t_end] ./
        (LP.LH_v0(earth_param_set) * 1000) .* (1e3 * 24 * 3600),
        color = :black,
    ) # not sure units

    # ax_SW_OUT
    p_SWOUT_m = lines!(
        ax_SWOUT,
        datetime2unix.(climaland.DateTime),
        climaland.SW_out,
        color = :blue,
    )
    p_SWOUT_d = lines!(
        ax_SWOUT,
        datetime2unix.(inputs.DateTime[index_t_start:index_t_end]),
        inputs.SW_OUT[index_t_start:index_t_end],
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
        ["Observations", "ClimaLand"],
        "",
        position = :rt,
        orientation = :horizontal,
    )

    ylims!(ax_C, (0, 40))
    ylims!(ax_W, (0, 30))
    ylims!(ax_SWOUT, (0, 200))

    [
        xlims!(
            axes,
            (
                datetime2unix(climaland.DateTime[1]),
                datetime2unix(climaland.DateTime[end]),
            ),
        ) for axes in [ax_C, ax_W, ax_SWOUT]
    ]

    fig
    return fig
end

# 2. Time series of SWC, Precip, moisture stress, stomatal conductance
function timeseries_H2O_fig(
    inputs,
    climaland,
    earth_param_set;
    dashboard = false,
) # will run for any inputs or climaland output of FLUXNET sites

    if dashboard == true
        WGLMakie.activate!() # for dashboards
    end

    index_t_start = findfirst(isequal(climaland.DateTime[1]), inputs.DateTime)
    index_t_end = findfirst(isequal(climaland.DateTime[end]), inputs.DateTime)


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

    # for time series, CairoMakie should allow DateTime type soon (but not yet)
    # so the 2 lines of code below are a trick to be able to use DateTime - will be removed later
    dateticks =
        optimize_ticks(climaland.DateTime[1], climaland.DateTime[end])[1][2:(end - 1)] # first and last are weirdly placed

    # add plots into axis ax_H2O
    p_H2O_m = lines!(
        ax_H2O,
        datetime2unix.(climaland.DateTime),
        climaland.θ_l,
        color = :green,
    )
    p_H2O_d = lines!(
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
    p_MS = lines!(
        ax_MS,
        datetime2unix.(climaland.DateTime),
        climaland.β,
        color = :green,
    ) # not sure about units

    # Stomatal conductance
    p_SC = lines!(
        ax_SC,
        datetime2unix.(climaland.DateTime),
        climaland.gs,
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
        ["Observations", "ClimaLand"],
        "",
        position = :rt,
        orientation = :horizontal,
    )

    #ylims!(ax_C, (0, 40))
    #ylims!(ax_W, (0, 30))
    #ylims!(ax_SWOUT, (0, 200))

    [
        xlims!(
            axes,
            (
                datetime2unix(climaland.DateTime[1]),
                datetime2unix(climaland.DateTime[end]),
            ),
        ) for axes in [ax_H2O, ax_H2O_rain, ax_MS, ax_SC]
    ]

    fig
    return fig
end

# 3. Fingerprint plot
function fingerprint_fig(inputs, climaland, earth_param_set; dashboard = false) # will run for any inputs or climaland output of FLUXNET sites

    if dashboard == true
        WGLMakie.activate!() # for dashboards
    end

    index_t_start = findfirst(isequal(climaland.DateTime[1]), inputs.DateTime)
    index_t_end = findfirst(isequal(climaland.DateTime[end]), inputs.DateTime)

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
    ax_W = Axis(
        fig[2, 1],
        ylabel = "Hour of the day",
        xlabel = "Date",
        title = L"\text{ET (mm)}",
    ) # h2o fluxes
    ax_M = Axis(
        fig[3, 1],
        ylabel = "Hour of the day",
        xlabel = "Date",
        title = L"\text{soil moisture}",
    ) # soil moisture
    ax_R = Axis(
        fig[4, 1],
        ylabel = "Hour of the day",
        xlabel = "Date",
        title = L"\text{incoming shortwave radiation}",
    ) # radiation

    # for time series, CairoMakie should allow DateTime type soon (but not yet)
    # so the 2 lines of code below are a trick to be able to use DateTime - will be removed later
    dateticks =
        optimize_ticks(inputs.DateTime[1], inputs.DateTime[end])[1][2:(end - 1)] # first and last are weirdly placed

    # Fingerprint plot
    hm_GPP = heatmap!(
        ax_C,
        datetime2unix.(DateTime.(Date.(inputs.DateTime))),
        hour.(inputs.DateTime) .+ (minute.(inputs.DateTime) ./ 60),
        inputs.GPP .* 1e6,
    )
    Colorbar(fig[1, 2], hm_GPP)

    hm_ET = heatmap!(
        ax_W,
        datetime2unix.(DateTime.(Date.(inputs.DateTime))),
        hour.(inputs.DateTime) .+ (minute.(inputs.DateTime) ./ 60),
        inputs.LE ./ (LP.LH_v0(earth_param_set) * 1000) .* (1e3 * 24 * 3600),
    )
    Colorbar(fig[2, 2], hm_ET)

    hm_M = heatmap!(
        ax_M,
        datetime2unix.(DateTime.(Date.(inputs.DateTime))),
        hour.(inputs.DateTime) .+ (minute.(inputs.DateTime) ./ 60),
        inputs.SWC,
    )
    Colorbar(fig[3, 2], hm_M)

    hm_R = heatmap!(
        ax_R,
        datetime2unix.(DateTime.(Date.(inputs.DateTime))),
        hour.(inputs.DateTime) .+ (minute.(inputs.DateTime) ./ 60),
        inputs.SW_IN,
    )
    Colorbar(fig[4, 2], hm_R)

    ax_C.xticks[] =
        (datetime2unix.(dateticks), Dates.format.(dateticks, "mm/dd"))

    ax_W.xticks[] =
        (datetime2unix.(dateticks), Dates.format.(dateticks, "mm/dd"))

    ax_M.xticks[] =
        (datetime2unix.(dateticks), Dates.format.(dateticks, "mm/dd"))

    ax_R.xticks[] =
        (datetime2unix.(dateticks), Dates.format.(dateticks, "mm/dd"))

    fig
    return fig
end

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
    diurnal_p = lines!(
        ax,
        0.5:1:23.5,
        getindex.(hourlyquantile[1:24], 2),
        color = color,
        linestyle = linestyle,
    )
    diurnal_q = band!(
        ax,
        0.5:1:23.5,
        getindex.(hourlyquantile[1:24], 1),
        getindex.(hourlyquantile[1:24], 3),
        color = (color, alpha),
    )
    return diurnal_p
end

function diurnals_fig(inputs, climaland, earth_param_set; dashboard = false) # will run for any inputs or climaland output of FLUXNET sites

    if dashboard == true
        WGLMakie.activate!() # for dashboards
    end

    index_t_start = findfirst(isequal(climaland.DateTime[1]), inputs.DateTime)
    index_t_end = findfirst(isequal(climaland.DateTime[end]), inputs.DateTime)


    fig = Figure(size = (1000, 1000))
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    ax_C = Axis(
        fig[1, 1],
        ylabel = L"\text{CO}_{2} \, (\mu\text{mol m}^{-2} \, \text{s}^{-1})",
    ) # C fluxes
    ax_W = Axis(fig[2, 1], ylabel = L"\text{H}_{2}\text{O} \, \text{(mm)}") # h2o fluxes
    ax_SIF = Axis(fig[3, 1], ylabel = L"\text{SIF}") # SIF
    ax_E = Axis(
        fig[4, 1],
        ylabel = L"\text{Radiation} \, (\text{W} \, \text{m}^{-2})",
        xlabel = L"\text{Hour of the day}",
        xgridvisible = false,
    ) # shortwave out

    # CO2 fluxes
    # model
    p_GPP_m = diurnal_plot!(
        fig,
        ax_C,
        climaland.DateTime,
        climaland.GPP .* 1e6,
        :green,
    )
    p_RA_m = diurnal_plot!(
        fig,
        ax_C,
        climaland.DateTime,
        climaland.Ra .* 1e6,
        :black,
    )
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
        climaland.DateTime,
        climaland.transpiration .* 1e3 .* 24 .* 3600,
        :blue,
    )
    # data
    p_ET_d = diurnal_plot!(
        fig,
        ax_W,
        inputs.DateTime[index_t_start:index_t_end],
        inputs.LE[index_t_start:index_t_end] ./
        (LP.LH_v0(earth_param_set) * 1000) .* (1e3 * 24 * 3600),
        :blue,
        alpha = 0.1,
        linestyle = :dot,
    )

    # SIF
    p_SIF_m =
        diurnal_plot!(fig, ax_SIF, climaland.DateTime, climaland.SIF, :black)


    # Energy fluxes
    # model
    # diurnal_plot!(fig, ax_E, climaland.DateTime, climaland.LW_out, :red)
    p_SWout_m =
        diurnal_plot!(fig, ax_E, climaland.DateTime, climaland.SW_out, :red)
    # data
    p_SWout_d = diurnal_plot!(
        fig,
        ax_E,
        inputs.DateTime[index_t_start:index_t_end],
        inputs.SW_OUT[index_t_start:index_t_end],
        :red,
        alpha = 0.1,
        linestyle = :dot,
    )

    [xlims!(axes, (0, 24)) for axes in [ax_C, ax_W, ax_SIF, ax_E]]

    axislegend(
        ax_C,
        [p_GPP_d, p_GPP_m, p_RA_m],
        ["GPP Obs.", "GPP model", "Ra model"],
        "",
        position = :rt,
        orientation = :vertical,
    )
    axislegend(
        ax_W,
        [p_ET_d, p_ET_m],
        ["ET obs.", "ET model"],
        "",
        position = :rt,
        orientation = :vertical,
    )
    axislegend(
        ax_SIF,
        [p_SIF_m],
        ["SIF model"],
        "",
        position = :rt,
        orientation = :vertical,
    )
    axislegend(
        ax_E,
        [p_SWout_d, p_SWout_m],
        ["SWout obs.", "SWout model"],
        "",
        position = :rt,
        orientation = :vertical,
    )

    hidexdecorations!(ax_C)
    hidexdecorations!(ax_W)
    hidexdecorations!(ax_SIF)

    fig
    return fig
end

# 5. Cumulative P and ET
function cumulative_H2O_fig(
    inputs,
    climaland,
    earth_param_set;
    dashboard = false,
)

    if dashboard == true
        WGLMakie.activate!() # for dashboards
    end

    index_t_start = findfirst(isequal(climaland.DateTime[1]), inputs.DateTime)
    index_t_end = findfirst(isequal(climaland.DateTime[end]), inputs.DateTime)


    fig = Figure(size = (1000, 1000))
    fontsize_theme = Theme(fontsize = 20)
    set_theme!(fontsize_theme)

    ax = Axis(fig[1, 1], ylabel = "Cumulative water flux (mm)")

    ET_m = climaland.transpiration .* 1e3 .* 24 .* 3600
    ET_obs =
        inputs.LE[index_t_start:index_t_end] ./
        (LP.LH_v0(earth_param_set) * 1000) .* (1e3 * 24 * 3600)
    P_obs = inputs.P[index_t_start:index_t_end] .* 1e3 .* 24 .* 3600

    p_P = lines!(ax, (1 / 48):(1 / 48):(length(P_obs) / 48), cumsum(P_obs)) # Precip
    p_ET_m = lines!(ax, (1 / 48):(1 / 48):(length(ET_m) / 48), cumsum(ET_m)) # ET model
    p_ET_d = lines!(ax, (1 / 48):(1 / 48):(length(ET_obs) / 48), cumsum(ET_obs)) # ET observations

    axislegend(
        ax,
        [p_P, p_ET_m, p_ET_d],
        ["Precipitation", "ET modeled", "ET observed"],
        "",
        position = :lt,
        orientation = :vertical,
    )

    fig
    return fig
end

# 6. Energy balance closure (L + H = Rn - G)

function make_plots(
    inputs,
    climaland;
    FT = Float64,
    save_fig = true,
    dashboard = false,
) # will run for any inputs or climaland output of FLUXNET sites

    if dashboard == true
        WGLMakie.activate!() # for dashboards
    else
        CairoMakie.activate!()
    end

    earth_param_set = LP.LandParameters(FT)
    args = (inputs, climaland, earth_param_set)

    fig1 = timeseries_fluxes_fig(args...)
    fig2 = timeseries_H2O_fig(args...)
    fig3 = fingerprint_fig(args...)
    fig4 = diurnals_fig(args...)
    fig5 = cumulative_H2O_fig(args...)

    if save_fig == true
        if isdir("figures")
            nothing
        else
            mkdir("figures")
        end

        names = [
            "timeseries_fluxes.pdf",
            "timeseries_H2O.pdf",
            "fingerprint.pdf",
            "diurnals.pdf",
            "cumulative_water.pdf",
        ]

        [
            save(joinpath("figures", name), fig) for
            (name, fig) in zip(names, [fig1, fig2, fig3, fig4, fig5])
        ]

        return nothing
    end

    if save_fig == false
        return (
            timeseries = fig1,
            water = fig2,
            fingerprint = fig3,
            diurnals = fig4,
            cumulative_water = fig5,
        )
    end

end

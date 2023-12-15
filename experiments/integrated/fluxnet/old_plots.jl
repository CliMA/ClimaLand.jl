
# Plotting
daily = sol.t ./ 3600 ./ 24
savedir = joinpath(climalsm_dir, "experiments/integrated/fluxnet/$site_ID/out/")

if !isdir(savedir)
    mkdir(savedir)
end

# Number of days to plot
num_days = N_days - N_spinup_days

# Plot model diurnal cycles without data comparisons
# Autotrophic Respiration
AR =
    [
        parent(sv.saveval[k].canopy.autotrophic_respiration.Ra)[1] for
        k in 1:length(sv.saveval)
    ] .* 1e6
plot_daily_avg("AutoResp", AR, dt * n, num_days, "μmol/m^2/s", savedir)

# Plot all comparisons of model diurnal cycles to data diurnal cycles
# GPP
model_GPP =
    [
        parent(sv.saveval[k].canopy.photosynthesis.GPP)[1] for
        k in 1:length(sv.saveval)
    ] .* 1e6

if drivers.GPP.status == absent
    plot_daily_avg("GPP", model_GPP, dt * n, num_days, "μmol/m^2/s", savedir)
else
    GPP_data =
        drivers.GPP.values[Int64(t_spinup ÷ DATA_DT):Int64(tf ÷ DATA_DT)] .* 1e6
    plot_avg_comp(
        "GPP",
        model_GPP,
        dt * n,
        GPP_data,
        FT(DATA_DT),
        num_days,
        drivers.GPP.units,
        savedir,
    )
end

# SW_OUT
SW_out_model = [parent(sv.saveval[k].SW_out)[1] for k in 1:length(sv.saveval)]
if drivers.SW_OUT.status == absent
    plot_daily_avg("SW OUT", SW_out_model, dt * n, num_days, "w/m^2", savedir)
else
    SW_out_data = FT.(drivers.SW_OUT.values)[Int64(t_spinup ÷ DATA_DT):Int64(
        tf ÷ DATA_DT,
    )]
    plot_avg_comp(
        "SW OUT",
        SW_out_model,
        dt * n,
        SW_out_data,
        FT(DATA_DT),
        num_days,
        drivers.SW_OUT.units,
        savedir,
    )
end

# LW_OUT
LW_out_model = [parent(sv.saveval[k].LW_out)[1] for k in 1:length(sv.saveval)]
if drivers.LW_OUT.status == absent
    plot_daily_avg("LW OUT", LW_out_model, dt * n, num_days, "w/m^2", savedir)
else
    LW_out_data = FT.(drivers.LW_OUT.values)[Int64(t_spinup ÷ DATA_DT):Int64(
        tf ÷ DATA_DT,
    )]
    plot_avg_comp(
        "LW OUT",
        LW_out_model,
        dt * n,
        LW_out_data,
        FT(DATA_DT),
        num_days,
        drivers.LW_OUT.units,
        savedir,
    )
end

# ET
T =
    [
        parent(sv.saveval[k].canopy.conductance.transpiration)[1] for
        k in 1:length(sol.t)
    ] .* (1e3 * 24 * 3600)
E =
    [parent(sv.saveval[k].soil_evap)[1] for k in 1:length(sol.t)] .* (1e3 * 24 * 3600)
ET_model = T .+ E
if drivers.LE.status == absent
    plot_daily_avg("ET", ET_model, dt * n, num_days, "mm/day", savedir)
else
    measured_T =
        drivers.LE.values ./ (LSMP.LH_v0(earth_param_set) * 1000) .*
        (1e3 * 24 * 3600)
    ET_data = measured_T[Int64(t_spinup ÷ DATA_DT):Int64(tf ÷ DATA_DT)]
    plot_avg_comp(
        "ET",
        ET_model,
        dt * n,
        ET_data,
        FT(DATA_DT),
        num_days,
        "mm/day",
        savedir,
    )
end

# Sensible Heat Flux
SHF_soil = [parent(sv.saveval[k].soil_shf)[1] for k in 1:length(sol.t)]
SHF_canopy =
    [parent(sv.saveval[k].canopy.energy.shf)[1] for k in 1:length(sol.t)]
SHF_model = SHF_soil + SHF_canopy
if drivers.H.status == absent
    plot_daily_avg("SHF", SHF_model, dt * n, num_days, "w/m^2", savedir)
else
    SHF_data = drivers.H.values[Int64(t_spinup ÷ DATA_DT):Int64(tf ÷ DATA_DT)]
    plot_avg_comp(
        "SHF",
        SHF_model,
        dt * n,
        SHF_data,
        FT(DATA_DT),
        N_days - N_spinup_days,
        drivers.H.units,
        savedir,
    )
end

# Ground Heat Flux
G_model = [
    (
        parent(sv.saveval[k].soil_shf)[1] + parent(sv.saveval[k].soil_lhf)[1] -
        parent(sv.saveval[k].soil_LW_n)[1] -
        parent(sv.saveval[k].soil_SW_n)[1]
    ) for k in 1:length(sol.t)
]
if drivers.G.status == absent
    plot_daily_avg("G", G_model, dt * n, num_days, "w/m^2", savedir)
else
    G_data =
        FT.(drivers.G.values)[Int64(t_spinup ÷ DATA_DT):Int64(tf ÷ DATA_DT)]
    plot_avg_comp(
        "G",
        G_model,
        dt * n,
        G_data,
        FT(DATA_DT),
        N_days - N_spinup_days,
        drivers.G.units,
        savedir,
    )
end

# Water stress factor
β = [parent(sv.saveval[k].canopy.hydraulics.β)[1] for k in 1:length(sol.t)]
plt1 =
    Plots.plot(size = (1500, 400), xlabel = "Day of year", margin = 10Plots.mm)
Plots.plot!(
    plt1,
    daily,
    β,
    label = "Model",
    xlim = [minimum(daily), maximum(daily)],
    title = "Moisture stress factor",
)
Plots.savefig(joinpath(savedir, "moisture_stress.png"))

# Stomatal conductance
g_stomata =
    [parent(sv.saveval[k].canopy.conductance.gs)[1] for k in 1:length(sol.t)]
plt1 =
    Plots.plot(size = (1500, 400), xlabel = "Day of year", margin = 10Plots.mm)
Plots.plot!(
    plt1,
    daily,
    g_stomata,
    label = "Model",
    xlim = [minimum(daily), maximum(daily)],
    title = "Stomatal conductance (mol/m^2/s)",
)
Plots.savefig(joinpath(savedir, "stomatal_conductance.png"))

# Soil water content
if drivers.SWC.status != absent
    # Current resolution has the first layer at 0.1 cm, the second at 5cm.
    plt1 = Plots.plot(size = (1500, 800))
    Plots.plot!(
        plt1,
        daily,
        [parent(sol.u[k].soil.ϑ_l)[end - 1] for k in 1:1:length(sol.t)],
        label = "5cm",
        xlim = [minimum(daily), maximum(daily)],
        ylim = [0.05, 0.55],
        xlabel = "Days",
        ylabel = "SWC [m/m]",
        color = "blue",
        margin = 10Plots.mm,
    )

    plot!(
        plt1,
        daily,
        [parent(sol.u[k].soil.θ_i)[end - 1] for k in 1:1:length(sol.t)],
        color = "cyan",
        label = "Ice, 5cm",
    )

    Plots.plot!(plt1, seconds ./ 3600 ./ 24, drivers.SWC.values, label = "Data")
    plt2 = Plots.plot(
        seconds ./ 3600 ./ 24,
        drivers.P.values .* (-1e3 * 24 * 3600),
        label = "Data",
        ylabel = "Precipitation [mm/day]",
        xlim = [minimum(daily), maximum(daily)],
        margin = 10Plots.mm,
        ylim = [-200, 0],
        size = (1500, 400),
    )
    Plots.plot(plt2, plt1, layout = grid(2, 1, heights = [0.2, 0.8]))
    Plots.savefig(joinpath(savedir, "soil_water_content.png"))
end

# Cumulative ET

dt_model = sol.t[2] - sol.t[1]
dt_data = seconds[2] - seconds[1]
# Find which index in the data our simulation starts at:
idx = argmin(abs.(seconds .- sol.t[1]))
if drivers.LE.status != absent
    Plots.plot(
        seconds ./ 24 ./ 3600,
        cumsum(measured_T[:]) * dt_data,
        label = "Data ET",
    )

    Plots.plot!(
        seconds ./ 24 ./ 3600,
        cumsum(drivers.P.values[:]) * dt_data * (1e3 * 24 * 3600),
        label = "Data P",
    )
    Plots.plot!(
        daily,
        cumsum(T .+ E) * dt_model .+ cumsum(measured_T[:])[idx] * dt_data,
        label = "Model ET",
    )

    Plots.plot!(
        ylabel = "∫ Water fluxes dt",
        xlabel = "Days",
        margins = 10Plots.mm,
    )
    Plots.savefig(joinpath(savedir, "cumul_p_et.png"))
end

# Soil Temperature
data_times = [0:DATA_DT:(num_days * S_PER_DAY);]
model_times = [0:(n * dt):(num_days * S_PER_DAY);]

# The second layer is ~ 5cm, third is at 11cm
soil_T_5 = [parent(sv.saveval[k].soil.T)[end - 1] for k in 1:length(sol.t)]
soil_T_5_avg = compute_diurnal_avg(soil_T_5, model_times, num_days)
soil_T_10 = [parent(sv.saveval[k].soil.T)[end - 2] for k in 1:length(sol.t)]
soil_T_10_avg = compute_diurnal_avg(soil_T_10, model_times, num_days)

TA_avg = compute_diurnal_avg(
    FT.(drivers.TA.values)[Int64(t_spinup ÷ DATA_DT):Int64(tf ÷ DATA_DT)],
    data_times,
    num_days,
)
if drivers.TS.status != absent
    TS_avg = compute_diurnal_avg(
        FT.(drivers.TS.values)[Int64(t_spinup ÷ DATA_DT):Int64(tf ÷ DATA_DT)],
        data_times,
        num_days,
    )
end

plt1 = Plots.plot(size = (1500, 400))
if drivers.TS.status != absent
    Plots.plot!(
        plt1,
        0.5:0.5:24,
        TS_avg,
        label = "Tsoil (data)",
        title = "Temperature",
    )
end
Plots.plot!(plt1, 0.5:0.5:24, TA_avg, label = "Tair (data)")
Plots.plot!(plt1, 0.5:0.5:24, soil_T_5_avg, label = "Tsoil (model; 5cm)")
Plots.plot!(plt1, 0.5:0.5:24, soil_T_10_avg, label = "Tsoil (model; 11cm)")
Plots.plot!(plt1, xlabel = "Hour of day", ylabel = "Average over Simulation")
Plots.plot!(plt1, margins = 10Plots.mm)
Plots.savefig(joinpath(savedir, "soil_temperature.png"))

# Temperatures
soil_T_sfc = [parent(sv.saveval[k].soil.T)[end] for k in 1:length(sol.t)]
soil_T_sfc_avg = compute_diurnal_avg(soil_T_sfc, model_times, num_days)

canopy_T = [
    parent(
        ClimaLSM.Canopy.canopy_temperature(
            land.canopy.energy,
            land.canopy,
            sol.u[k],
            sv.saveval[k],
            sol.t[k],
        ),
    )[1] for k in 1:length(sol.t)
]
canopy_T_avg = compute_diurnal_avg(canopy_T, model_times, num_days)

plt1 = Plots.plot(size = (1500, 400))
if drivers.TS.status != absent
    Plots.plot!(
        plt1,
        0.5:0.5:24,
        TS_avg,
        label = "Soil-D",
        title = "Temperature",
    )
end
Plots.plot!(plt1, 0.5:0.5:24, TA_avg, label = "Atmos-D")

Plots.plot!(plt1, 0.5:0.5:24, soil_T_sfc_avg, label = "Soil-M-2.5cm")

Plots.plot!(plt1, 0.5:0.5:24, canopy_T_avg, label = "Canopy-M")
Plots.plot!(plt1, xlabel = "Hour of day", ylabel = "Average over Simulation")
Plots.plot!(plt1, margins = 10Plots.mm)
Plots.savefig(joinpath(savedir, "temperature.png"))

# Run script with comand line argument "save" to save model output to CSV
if length(ARGS) ≥ 1 && ARGS[1] == "save"
    # Formats fields as semicolon seperated strings
    field_to_array = (field) -> join(parent(field), ';')
    # Recursively unpacks a nested NamedTuple of fields into an array of strings
    function unpack(tup, data)
        for entry in tup
            if entry isa NamedTuple
                unpack(entry, data)
            else
                push!(data, field_to_array(entry))
            end
        end
    end
    # Recursively extracts the names of all fields in a nested namedTuple
    function extract_names(nt, names)
        for entry in pairs(nt)
            if entry[2] isa NamedTuple
                extract_names(entry[2], names)
            else
                push!(names, entry[1])
            end
        end
    end
    # Collect unpacked data from each timestep into an array
    timestamps = [[]]
    push!(timestamps[1], "Timestep")
    extract_names(sv.saveval[1], timestamps[1])
    local cnt = 0
    for timestamp in sv.saveval
        cnt = cnt + 1
        save_data = Any[cnt]
        unpack(timestamp, save_data)
        push!(timestamps, save_data)
    end
    # Write all data to a csv file
    writedlm(joinpath(savedir, "model_output.csv"), timestamps, ',')
    @info "Saved model output to $(savedir)model_output.csv"
end

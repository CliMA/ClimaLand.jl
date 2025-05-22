import SciMLBase
import ClimaComms
ClimaComms.@import_required_backends
import ClimaTimeSteppers as CTS
import ClimaCore
@show pkgversion(ClimaCore)
using ClimaUtilities.ClimaArtifacts

import ClimaDiagnostics
import ClimaAnalysis
import ClimaAnalysis.Visualize as viz

import ClimaUtilities
import ClimaUtilities.OnlineLogging: WallTimeInfo, report_walltime
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.TimeManager: ITime
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP

using ClimaLand
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP

using Statistics
using GeoMakie
using CairoMakie
using Dates
import NCDatasets

using Poppler_jll: pdfunite
greet = true
const FT = Float64;
time_interpolation_method = LinearInterpolation(PeriodicCalendar())
context = ClimaComms.context()
ClimaComms.init(context)
device = ClimaComms.device()
device_suffix = device isa ClimaComms.CPUSingleThreaded ? "cpu" : "gpu"
root_path = "california_longrun_$(device_suffix)"
diagnostics_outdir = joinpath(root_path, "regional_diagnostics")
outdir =
    ClimaUtilities.OutputPathGenerator.generate_output_path(diagnostics_outdir)

    start_date = DateTime(2008)
    t0 = 0.0
    tf = 60 * 60.0 * 24 * 340
    Δt =900.0
    nelements = (5, 5, 15)
    if greet
        @info "Run: Regional Soil-Canopy-Snow Model"
        @info "Resolution: $nelements"
        @info "Timestep: $Δt s"
        @info "Duration: $(tf - t0) s"
    end

    t0 = ITime(t0, epoch = start_date)
    tf = ITime(tf, epoch = start_date)
    Δt = ITime(Δt, epoch = start_date)
    t0, tf, Δt = promote(t0, tf, Δt)
    earth_param_set = LP.LandParameters(FT)
    radius = FT(6378.1e3)
    depth = FT(10)
    center_long, center_lat = FT(-76), FT(3)
    delta_m = FT(100_000) # in meters, this is about a 2 degree simulation
    domain = ClimaLand.Domains.HybridBox(;
        xlim = (delta_m, delta_m),
        ylim = (delta_m, delta_m),
        zlim = (-depth, FT(0)),
        nelements = nelements,
        longlat = (center_long, center_lat),
        dz_tuple = FT.((2.0, 0.05)),
    );
    surface_space = domain.space.surface;
    subsurface_space = domain.space.subsurface;

    # Forcing data
    era5_ncdata_path =
        ClimaLand.Artifacts.era5_land_forcing_data2008_path(; context)
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        era5_ncdata_path,
        surface_space,
        start_date,
        earth_param_set,
        FT,
    );
    spatially_varying_soil_params =
        ClimaLand.default_spatially_varying_soil_parameters(
            subsurface_space,
            surface_space,
            FT,
        )
    (;
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm,
        K_sat,
        S_s,
        θ_r,
        PAR_albedo_dry,
        NIR_albedo_dry,
        PAR_albedo_wet,
        NIR_albedo_wet,
        f_max,
    ) = spatially_varying_soil_params
    soil_params = Soil.EnergyHydrologyParameters(
        FT;
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm,
        K_sat,
        S_s,
        θ_r,
        PAR_albedo_dry = PAR_albedo_dry,
        NIR_albedo_dry = NIR_albedo_dry,
        PAR_albedo_wet = PAR_albedo_wet,
        NIR_albedo_wet = NIR_albedo_wet,
    );

    f_over = FT(3.28) # 1/m
    R_sb = FT(1.484e-4 / 1000) # m/s
    runoff_model = ClimaLand.Soil.Runoff.TOPMODELRunoff{FT}(;
        f_over = f_over,
        f_max = f_max,
        R_sb = R_sb,
    );

    soil_args = (domain = domain, parameters = soil_params)
    soil_model_type = Soil.EnergyHydrology{FT}
    sources = (Soil.PhaseChange{FT}(),)# sublimation and subsurface runoff are added automatically
    top_bc = ClimaLand.Soil.AtmosDrivenFluxBC(
        atmos,
        radiation,
        runoff_model,
        (:soil,),
    )
    zero_flux = Soil.HeatFluxBC((p, t) -> 0.0)
    boundary_conditions = (;
        top = top_bc,
        bottom = Soil.EnergyFreeDrainage()
	)
    soil = soil_model_type(;
        boundary_conditions = boundary_conditions,
        sources = sources,
        soil_args...,
    );

    Y, p, cds = initialize(soil);

    ic_path = ClimaLand.Artifacts.soil_ic_2008_50m_path(; context = context)
    evaluate!(p.drivers.T, atmos.T, t0)
    T_bounds = extrema(p.drivers.T)

    ClimaLand.Simulations.set_soil_initial_conditions!(
        Y,
        ν,
        θ_r,
        subsurface_space,
        ic_path,
        soil,
        T_bounds,
    )
    set_initial_cache! = make_set_initial_cache(soil)
    exp_tendency! = make_exp_tendency(soil)
    imp_tendency! = ClimaLand.make_imp_tendency(soil)
    jacobian! = ClimaLand.make_jacobian(soil)
    set_initial_cache!(p, Y, t0)

    # set up jacobian info
    jac_kwargs = (;
        jac_prototype = ClimaLand.FieldMatrixWithSolver(Y),
        Wfact = jacobian!,
    )

    prob = SciMLBase.ODEProblem(
        CTS.ClimaODEFunction(
            T_exp! = exp_tendency!,
            T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
            dss! = ClimaLand.dss!,
        ),
        Y,
        (t0, tf),
        p,
    );

    updateat = [promote(t0:(ITime(3600 * 3)):tf...)...];
    drivers = ClimaLand.get_drivers(soil);
    updatefunc = ClimaLand.make_update_drivers(drivers);

    # ClimaDiagnostics

    t_first_save = ITime(0, epoch = start_date)
    saveat = [promote(t_first_save:(ITime(900.0)):tf...)...]
    sv = (;
        t = Array{ITime{Int64, Second, DateTime}}(undef, length(saveat)),
        saveval = Array{NamedTuple}(undef, length(saveat)),
    )
    saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat)
    
    driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc);

    nancheck_freq = Dates.Month(1)
    nancheck_cb = ClimaLand.NaNCheckCallback(nancheck_freq, start_date, Δt)
    report_cb = ClimaLand.ReportCallback(1000)
    cb = SciMLBase.CallbackSet(driver_cb, saving_cb, nancheck_cb, report_cb);


    # Define timestepper and ODE algorithm
    stepper = CTS.ARS111()
    ode_algo = CTS.IMEXAlgorithm(
        stepper,
        CTS.NewtonsMethod(
            max_iters = 3,
            update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        ),
    )
    sol = SciMLBase.solve(prob, ode_algo; dt = Δt, callback = cb, saveat = saveat, adaptive = false);
  

extract2d(x) = mean(Array(parent(x))[1, 1, 1, 9])
extract3d(x; i = 1) =
    Array(parent(x))[i,1, 1, 1, 9]
extract_column(x) = Array(parent(x))[:,1, 1, 1, 9]
id_last = 29278 + 1000
id_first = 1
S = [
    extract3d((sv.saveval[k].soil.θ_l .- θ_r) ./ (ν .- θ_r)) for k in id_first:id_last
];
ϑ_l = [extract3d(sol.u[k].soil.ϑ_l) for k in id_first:id_last];
θ_l = [extract3d(sv.saveval[k].soil.θ_l) for k in id_first:id_last];
ν_eff = [extract3d(ν .- sol.u[k].soil.θ_i) for k in id_first:id_last];

ψ = [extract3d(sv.saveval[k].soil.ψ) for k in id_first:id_last];
T_soil = [extract3d(sv.saveval[k].soil.T) for k in id_first:id_last];
T_soil2 = [extract3d(sv.saveval[k].soil.T; i = 2) for k in id_first:id_last];
T_soil3 = [extract3d(sv.saveval[k].soil.T; i = 3) for k in id_first:id_last];
T_soil4 = [extract3d(sv.saveval[k].soil.T; i = 4) for k in id_first:id_last];
T_soil5 = [extract3d(sv.saveval[k].soil.T; i = 5) for k in id_first:id_last];
T_soil6 = [extract3d(sv.saveval[k].soil.T; i = 6) for k in id_first:id_last];


θ_i = [extract3d(sol.u[k].soil.θ_i) for k in id_first:id_last];
θ_i2 = [extract3d(sol.u[k].soil.θ_i; i = 2) for k in id_first:id_last];
θ_i3 = [extract3d(sol.u[k].soil.θ_i; i = 3) for k in id_first:id_last];
θ_i4 = [extract3d(sol.u[k].soil.θ_i; i = 4) for k in id_first:id_last];
θ_i5 = [extract3d(sol.u[k].soil.θ_i; i = 5) for k in id_first:id_last];
θ_i6 = [extract3d(sol.u[k].soil.θ_i; i = 6) for k in id_first:id_last];

e1 = [extract3d(sol.u[k].soil.ρe_int) for k in id_first:id_last];
e2 = [extract3d(sol.u[k].soil.ρe_int; i = 2) for k in id_first:id_last];
e3 = [extract3d(sol.u[k].soil.ρe_int; i = 3) for k in id_first:id_last];
e4 = [extract3d(sol.u[k].soil.ρe_int; i = 4) for k in id_first:id_last];
e5 = [extract3d(sol.u[k].soil.ρe_int; i = 5) for k in id_first:id_last];
e6 = [extract3d(sol.u[k].soil.ρe_int; i = 6) for k in id_first:id_last];
top_bc = [extract2d(sv.saveval[k].soil.top_bc.water) for k in id_first:id_last]
top_bc_heat = [extract2d(sv.saveval[k].soil.top_bc.heat) for k in id_first:id_last]

precip = [extract2d(sv.saveval[k].drivers.P_liq) for k in id_first:id_last]
evap = [
    extract2d(sv.saveval[k].soil.turbulent_fluxes.vapor_flux_liq) for
    k in id_first:id_last
]
    lhf = [
    extract2d(sv.saveval[k].soil.turbulent_fluxes.lhf) for
    k in id_first:id_last
]
    shf = [
    extract2d(sv.saveval[k].soil.turbulent_fluxes.shf) for
    k in id_first:id_last
]
    rn = [extract2d(sv.saveval[k].soil.R_n) for k in id_first:id_last];

sub = [
    extract2d(sv.saveval[k].soil.turbulent_fluxes.vapor_flux_ice) for
    k in id_first:id_last
]
runoff = [
    extract2d(sv.saveval[k].soil.R_s) for
    k in id_first:id_last
]
runoff_ss = [
    extract2d(sv.saveval[k].soil.R_ss) for
    k in id_first:id_last
]

T_air = [extract2d(sv.saveval[k].drivers.T) for k in id_first:id_last]

days = [sv.t[k].counter for k in id_first:id_last] ./ 24 ./ 3600

S14 = [
           extract3d((sv.saveval[k].soil.θ_l .- θ_r) ./ (ν .- θ_r); i = 14) for k in id_first:id_last
       ];
S13 = [
           extract3d((sv.saveval[k].soil.θ_l .- θ_r) ./ (ν .- θ_r); i = 13) for k in id_first:id_last
       ];
S12 = [
           extract3d((sv.saveval[k].soil.θ_l .- θ_r) ./ (ν .- θ_r); i = 12) for k in id_first:id_last
       ];

S11 = [
           extract3d((sv.saveval[k].soil.θ_l .- θ_r) ./ (ν .- θ_r); i = 11) for k in id_first:id_last
       ];

S10 = [
           extract3d((sv.saveval[k].soil.θ_l .- θ_r) ./ (ν .- θ_r); i = 10) for k in id_first:id_last
       ];

S9 = [
           extract3d((sv.saveval[k].soil.θ_l .- θ_r) ./ (ν .- θ_r); i = 9) for k in id_first:id_last
       ];

S8 = [
           extract3d((sv.saveval[k].soil.θ_l .- θ_r) ./ (ν .- θ_r); i = 8) for k in id_first:id_last
       ];

S7 = [
           extract3d((sv.saveval[k].soil.θ_l .- θ_r) ./ (ν .- θ_r); i = 7) for k in id_first:id_last
       ];
       S6 = [
           extract3d((sv.saveval[k].soil.θ_l .- θ_r) ./ (ν .- θ_r); i = 6) for k in id_first:id_last
       ];
S5 = [
           extract3d((sv.saveval[k].soil.θ_l .- θ_r) ./ (ν .- θ_r); i = 5) for k in id_first:id_last
       ];
       S4 = [
           extract3d((sv.saveval[k].soil.θ_l .- θ_r) ./ (ν .- θ_r); i = 4) for k in id_first:id_last
       ];
S3 = [
           extract3d((sv.saveval[k].soil.θ_l .- θ_r) ./ (ν .- θ_r); i = 3) for k in id_first:id_last
       ];
       S2 = [
           extract3d((sv.saveval[k].soil.θ_l .- θ_r) ./ (ν .- θ_r); i = 2) for k in id_first:id_last
       ];
S1 = [
           extract3d((sv.saveval[k].soil.θ_l .- θ_r) ./ (ν .- θ_r); i = 1) for k in id_first:id_last
       ];

fig = CairoMakie.Figure(size = (1500, 1500))
ax = CairoMakie.Axis(fig[1, 1], xlabel = "Time", ylabel = "Water Content")
lines!(ax, days, ϑ_l .- extract3d(θ_r), label = "ϑ_l- θ_r")
lines!(ax, days, θ_i, label = "θ_i, 1")
lines!(ax, days, θ_i2, label = "θ_i, 2")
lines!(ax, days, θ_i3, label = "θ_i, 3")
lines!(ax, days, θ_i4, label = "θ_i, 4")
lines!(ax, days, θ_i5, label = "θ_i, 5")
lines!(ax, days, θ_i6, label = "θ_i, 6")

lines!(ax, days, S, label = "S_l")
axislegend(ax; position = :rt)

ax = CairoMakie.Axis(
    fig[2, 2],
    xlabel = "Time",
    ylabel = "Vol internal energy")
#    ylabel = "Matric Potential",
#    yscale = log10,
#)
#lines!(ax, days, abs.(ψ), label = "ψ(m)")
lines!(ax, days, e1, label = "e1")
lines!(ax, days, e2, label = "e2")
lines!(ax, days, e3, label = "e3")
lines!(ax, days, e4, label = "e4")
lines!(ax, days, e5, label = "e5")
lines!(ax, days, e6, label = "e6")

axislegend(ax)
ax = CairoMakie.Axis(fig[2, 3], xlabel = "Time", ylabel = "Temperature")
lines!(ax, days, T_soil, label = "T soil, 1")
lines!(ax, days, T_soil2, label = "T soil, 2")
lines!(ax, days, T_soil3, label = "T soil, 3")
lines!(ax, days, T_soil4, label = "T soil, 4")
lines!(ax, days, T_soil5, label = "T soil, 5")
lines!(ax, days, T_soil6, label = "T soil, 6")
lines!(ax, days, T_air[id_first:id_last], label = "T air")
axislegend(ax)
ax = CairoMakie.Axis(fig[2, 1], xlabel = "Time", ylabel = "Water Fluxes")
lines!(ax, days, evap, label = "E")
lines!(ax, days, sub, label = "S")
lines!(ax, days, runoff, label = "R_s")
lines!(ax, days, precip[id_first:id_last], label = "precip")
lines!(ax, days, top_bc, label = "BC")
axislegend(ax; position = :lb)

ax = CairoMakie.Axis(fig[1, 2], xlabel = "Time", ylabel = "Heat Fluxes")
lines!(ax, days, lhf, label = "lhf")
lines!(ax, days, shf, label = "shf")
lines!(ax, days, rn, label = "rn")
lines!(ax, days, top_bc_heat, label = "BC")
axislegend(ax)
ax = CairoMakie.Axis(fig[1, 3], xlabel = "Time", ylabel = "subsurface")
lines!(ax, days, S1, label = "S1")
lines!(ax, days, S2, label = "S2")
lines!(ax, days, S3, label = "S3")
lines!(ax, days, S4, label = "S4")
lines!(ax, days, S5, label = "S5")
lines!(ax, days, S6, label = "S6")

axislegend(ax; position=:lb)
CairoMakie.save("soil_sa_zoom_bot_bc.png", fig)

# # Regional run of soil model

import SciMLBase
import ClimaComms
ClimaComms.@import_required_backends
import ClimaTimeSteppers as CTS
using ClimaCore
using ClimaUtilities.ClimaArtifacts
import Interpolations
using Insolation

import ClimaDiagnostics
import ClimaAnalysis
import ClimaAnalysis.Visualize as viz
import ClimaUtilities


import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaParams as CP

using ClimaLand
using ClimaLand.Soil
using ClimaLand.Snow
using ClimaLand.Canopy
import ClimaLand
import ClimaLand.Parameters as LP

using Statistics
using CairoMakie
import GeoMakie
using Dates
import NCDatasets

const FT = Float64;

time_interpolation_method = LinearInterpolation(PeriodicCalendar())
context = ClimaComms.context()
device = ClimaComms.device()
t0 = 0.0
tf = 60 * 60.0 * 24 * 20
Δt = 450.0
nelements = (10, 10, 15)

earth_param_set = LP.LandParameters(FT)
depth = FT(50)
center_long, center_lat = FT(36), FT(-8)
delta_m = FT(200_000) #~ few x few degrees
domain = ClimaLand.Domains.HybridBox(;
    xlim = (delta_m, delta_m),
    ylim = (delta_m, delta_m),
    zlim = (-depth, FT(0)),
    nelements = nelements,
    npolynomial = 1,
    longlat = (center_long, center_lat),
    dz_tuple = FT.((10.0, 0.05)),
)
surface_space = domain.space.surface
subsurface_space = domain.space.subsurface

start_date = DateTime(2008)
# Forcing data
era5_artifact_path =
    ClimaLand.Artifacts.era5_land_forcing_data2008_folder_path(; context)
era5_ncdata_path = joinpath(era5_artifact_path, "era5_2008_1.0x1.0.nc")
atmos, radiation = ClimaLand.prescribed_forcing_era5(
    era5_ncdata_path,
    surface_space,
    start_date,
    earth_param_set,
    FT;
    time_interpolation_method = time_interpolation_method,
)
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
)
f_over = FT(3.28) # 1/m
R_sb = FT(1.484e-4 / 1000) # m/s
runoff_model = ClimaLand.Soil.Runoff.TOPMODELRunoff{FT}(;
    f_over = f_over,
    f_max = f_max,
    R_sb = R_sb,
)

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
    bottom = Soil.WaterHeatBC(;
        water = Soil.FreeDrainage(),
        heat = zero_flux,
    ),
)
soil = soil_model_type(;
    boundary_conditions = boundary_conditions,
    sources = sources,
    soil_args...,
)

Y, p, cds = initialize(soil)
z = ClimaCore.Fields.coordinate_field(cds.subsurface).z
lat = ClimaCore.Fields.coordinate_field(cds.subsurface).lat
# This function approximates the hydrostatic equilibrium solution in
# the vadose and unsaturated regimes by solving for ∂(ψ+z)/∂z = 0,
# assuming the transition between the two is at a coordinate of z_∇.

# The approximation arises because the porosity, residual water fraction,
# and van Genuchtem parameters are spatially varying but treated constant
# in solving for equilibrium. Finally, we make a plausible but total guess
# for the water table depth using the TOPMODEL parameters.
function hydrostatic_profile(
    lat::FT,
    z::FT,
    ν::FT,
    θ_r::FT,
    α::FT,
    n::FT,
    S_s::FT,
    fmax,
) where {FT}
    m = 1 - 1 / n
    zmin = FT(-50.0)
    zmax = FT(0.0)

    z_∇ = FT(zmin / 5.0 + (zmax - zmin) / 2.5 * (fmax - 0.35) / 0.7)
    if z > z_∇
        S = FT((FT(1) + (α * (z - z_∇))^n)^(-m))
        ϑ_l = S * (ν - θ_r) + θ_r
    else
        ϑ_l = -S_s * (z - z_∇) + ν
    end
    return FT(ϑ_l)
end
vg_α = hydrology_cm.α
vg_n = hydrology_cm.n
Y.soil.ϑ_l .= hydrostatic_profile.(lat, z, ν, θ_r, vg_α, vg_n, S_s, f_max)
Y.soil.θ_i .= FT(0.0)
T = FT(276.85)
ρc_s =
    Soil.volumetric_heat_capacity.(
        Y.soil.ϑ_l,
        Y.soil.θ_i,
        soil_params.ρc_ds,
        soil_params.earth_param_set,
    )
Y.soil.ρe_int .=
    Soil.volumetric_internal_energy.(
        Y.soil.θ_i,
        ρc_s,
        T,
        soil_params.earth_param_set,
    )

set_initial_cache! = make_set_initial_cache(soil)
exp_tendency! = make_exp_tendency(soil)
imp_tendency! = ClimaLand.make_imp_tendency(soil)
jacobian! = ClimaLand.make_jacobian(soil)
set_initial_cache!(p, Y, t0)

# set up jacobian info
jac_kwargs =
    (; jac_prototype = ImplicitEquationJacobian(Y), Wfact = jacobian!)

prob = SciMLBase.ODEProblem(
    CTS.ClimaODEFunction(
        T_exp! = exp_tendency!,
        T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
        dss! = ClimaLand.dss!,
    ),
    Y,
    (t0, tf),
    p,
)

updateat = Array(t0:(3600 * 3):tf)
drivers = ClimaLand.get_drivers(soil)
updatefunc = ClimaLand.make_update_drivers(drivers)
t_first_save = t0 #+ 395*24*3600.0
dt_save = Δt
saveat = Array(t_first_save:dt_save:tf)
sv = (;
    t = Array{Float64}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
)
saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat)

driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
cb = SciMLBase.CallbackSet(driver_cb, saving_cb)


# Define timestepper and ODE algorithm
stepper = CTS.ARS111()
ode_algo = CTS.IMEXAlgorithm(
    stepper,
    CTS.NewtonsMethod(
        max_iters = 6,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
    ),
)
sol = SciMLBase.solve(
    prob,
    ode_algo;
    dt = Δt,
    callback = cb,
    adaptive = false,
    saveat = saveat,
);

# found by looking at different inds of isnan.(Array(parent(ClimaLand.top_center_to_surface(sol.u[end].soil.θ_i))))[:,:,1,3]
extract2d(x) = mean(Array(parent(x))[1,1,1,3])
extract3d(x) =
    mean(Array(parent(ClimaLand.top_center_to_surface(x)))[1, 1, 1, 3])
# id_last = length(sol.t)
id_last = Int64((60 * 60 * 24 * 20) / Δt)
S = [
    extract3d((sv.saveval[k].soil.θ_l .- θ_r) ./ (ν .- θ_r)) for k in 1:id_last
];
ϑ_l = [extract3d(sol.u[k].soil.ϑ_l) for k in 1:id_last];

# Δz_soil  = [extract2d(sv.saveval[k].effective_soil_sfc_depth) for k in 1:id_last];
# T_eff_soil  = [extract2d(sv.saveval[k].effective_soil_sfc_T) for k in 1:id_last];

# ghf = [extract2d(sv.saveval[k].ground_heat_flux) for k in 1:id_last];
κ_soil = [extract3d(sv.saveval[k].soil.κ  ) for k in 1:id_last];
ψ = [extract3d(sv.saveval[k].soil.ψ) for k in 1:id_last];
T_soil = [extract3d(sv.saveval[k].soil.T) for k in 1:id_last];

θ_i = [extract3d(sol.u[k].soil.θ_i) for k in 1:id_last];
top_bc = [extract2d(sv.saveval[k].soil.top_bc.water) for k in 1:id_last]
top_bc_heat = [extract2d(sv.saveval[k].soil.top_bc.heat) for k in 1:id_last]

precip = [extract2d(sv.saveval[k].drivers.P_liq) for k in 1:id_last]
evap = [
    extract2d(sv.saveval[k].soil.turbulent_fluxes.vapor_flux_liq) for
    k in 1:id_last
]
sub = [
    extract2d(sv.saveval[k].soil.turbulent_fluxes.vapor_flux_ice) for
    k in 1:id_last
]

T_air = [extract2d(sv.saveval[k].drivers.T) for k in 1:id_last]

days = sv.t[1:id_last] ./ 24 ./ 3600


output_dir = joinpath(pwd(), "soil_region_savingcb_output_12-2")
!isdir(output_dir) && mkdir(output_dir)

using DelimitedFiles

# Get all state and cache variables we accessed above
vardict = ("S" => S,
    "ϑ_l" => ϑ_l,
    # "Δz_soil" => Δz_soil,
    # "T_eff_soil" => T_eff_soil,
    # "ghf" => ghf,
    "κ_soil" => κ_soil,
    "ψ" => ψ,
    "T_soil" => T_soil,
    "θ_i" => θ_i,
    "top_bc" => top_bc,
    "top_bc_heat" => top_bc_heat,
    "precip" => precip,
    "evap" => evap,
    "sub" => sub,
    "T_air" => T_air)

# Save everything to comma-delimited text files
for (varname, var) in vardict
    open(
        joinpath(output_dir, "$varname.txt"),
        "w",
    ) do io
        writedlm(io, var, ',')
    end;
end

fig = CairoMakie.Figure(size = (1500, 1500))
ax = CairoMakie.Axis(fig[1, 1], xlabel = "Time", ylabel = "Water Content")
#lines!(ax, days, ϑ_l .- extract3d(θ_r), label = "ϑ_l- θ_r")
lines!(ax, days, θ_i, label = "θ_i")
lines!(ax, days, S, label = "S_l")
axislegend(ax, position=:lb)



ax = CairoMakie.Axis(
    fig[2, 2],
    xlabel = "Time",
    ylabel = "Matric Potential",
    yscale = log10,
)
lines!(ax, days, abs.(ψ), label = "ψ(m)")
axislegend(ax, position=:lb)
ax = CairoMakie.Axis(fig[2, 3], xlabel = "Time", ylabel = "Temperature")
lines!(ax, days, T_soil, label = "T soil")
# lines!(ax, days, T_eff_soil, label = "T eff_soil")
lines!(ax, days, T_air, label = "T air")
axislegend(ax, position=:lb)
ax = CairoMakie.Axis(fig[2, 1], xlabel = "Time", ylabel = "Water Fluxes")
lines!(ax, days, evap, label = "E")
lines!(ax, days, sub, label = "S")
lines!(ax, days, precip, label = "precip")
lines!(ax, days, top_bc, label = "BC")
axislegend(ax, position=:lb)

ax = CairoMakie.Axis(fig[1, 2], xlabel = "Time", ylabel = "Heat Fluxes")
# lines!(ax, days, ghf, label = "ghf")
lines!(ax, days, top_bc_heat, label = "BC")
axislegend(ax, position=:lb)
CairoMakie.save(joinpath(output_dir, "soil.png"), fig)

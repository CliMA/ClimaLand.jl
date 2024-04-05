using CairoMakie
using Statistics
using ArtifactWrappers
using Dates
import SciMLBase
import ClimaTimeSteppers as CTS
using ClimaCore
using ClimaUtilities.ClimaArtifacts
import Interpolations

import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput

import ClimaUtilities.DataHandling
import NCDatasets
import ClimaParams as CP
using ClimaComms

using ClimaLand
using ClimaLand.Soil
import ClimaLand
import ClimaLand.Parameters as LP


context = ClimaComms.context()
outdir = joinpath(pkgdir(ClimaLand), "experiments/standalone/Soil/artifacts")
!ispath(outdir) && mkpath(outdir)
FT = Float64
radius = FT(6378.1e3);
depth = FT(50)
domain = ClimaLand.Domains.SphericalShell(;
    radius = radius,
    depth = depth,
    nelements = (101, 15),
    npolynomial = 1,
    dz_tuple = FT.((10.0, 0.1)),
);
surface_space = domain.space.surface
subsurface_space = domain.space.subsurface
# Read in f_max data and land sea mask
topmodel_dataset = ArtifactWrapper(
    @__DIR__,
    "processed_topographic_index 2.5 deg",
    ArtifactFile[ArtifactFile(
        url = "https://caltech.box.com/shared/static/dwa7g0uzhxd50a2z3thbx3a12n0r887s.nc",
        filename = "means_2.5_new.nc",
    ),],
)
infile_path = joinpath(get_data_folder(topmodel_dataset), "means_2.5_new.nc")
outfile_root =
    joinpath(pkgdir(ClimaLand), "experiments/standalone/Soil/static_data_cgll")

f_max = SpaceVaryingInput(
    infile_path,
    "fmax",
    surface_space;
    regridder_type = :InterpolationsRegridder,
)
mask = SpaceVaryingInput(
    infile_path,
    "landsea_mask",
    surface_space;
    regridder_type = :InterpolationsRegridder,
)

oceans_to_zero(field, mask) = mask > 0.5 ? field : eltype(field)(0)
f_over = FT(3.28) # 1/m
R_sb = FT(1.484e-4 / 1000) # m/s
runoff_model = ClimaLand.Soil.Runoff.TOPMODELRunoff{FT}(;
    f_over = f_over,
    f_max = f_max,
    R_sb = R_sb,
)
vg_α = ClimaCore.Fields.ones(subsurface_space) .* FT(0.2)
hydrology_cm = map(vg_α) do (α)
    FT = typeof(α)
    ClimaLand.Soil.vanGenuchten{FT}(; α = α, n = α + FT(2))
end
θ_r = ClimaCore.Fields.zeros(subsurface_space)
ν = ClimaCore.Fields.zeros(subsurface_space) .+ FT(0.5)
K_sat = ClimaCore.Fields.zeros(subsurface_space) .+ FT(1e-6)
S_s = ClimaCore.Fields.zeros(subsurface_space) .+ FT(1e-3)

soil_params = ClimaLand.Soil.RichardsParameters(;
    hydrology_cm = hydrology_cm,
    ν = ν,
    K_sat = K_sat,
    S_s = S_s,
    θ_r = θ_r,
)

era5_artifact_path = "/groups/esm/ClimaArtifacts/artifacts/era5_land_forcing_data2021"# @clima_artifact("era5_land_forcing_data2021")
ref_time = DateTime(2021);
t_start = 0.0
# Precipitation:
precip = TimeVaryingInput(
    joinpath(era5_artifact_path, "era5_2021_0.9x1.25.nc"),
    "tp",
    surface_space;
    reference_date = ref_time,
    t_start,
    regridder_type = :InterpolationsRegridder,
)
atmos = ClimaLand.PrescribedPrecipitation{FT, typeof(precip)}(precip)
bottom_bc = ClimaLand.Soil.WaterFluxBC((p, t) -> 0.0)
bc = (;
    top = ClimaLand.Soil.RichardsAtmosDrivenFluxBC(atmos, runoff_model),
    bottom = bottom_bc,
)
model = ClimaLand.Soil.RichardsModel{FT}(;
    parameters = soil_params,
    domain = domain,
    boundary_conditions = bc,
    sources = (),
    lateral_flow = false,
)
Y, p, t = initialize(model)
z = ClimaCore.Fields.coordinate_field(domain.space.subsurface).z
lat = ClimaCore.Fields.coordinate_field(domain.space.subsurface).lat
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
t0 = 0.0
tf = 3600.0 * 24
dt = 1800.0
Y.soil.ϑ_l .= hydrostatic_profile.(lat, z, ν, θ_r, vg_α, vg_α .+ 2, S_s, f_max)
@. Y.soil.ϑ_l = oceans_to_zero(Y.soil.ϑ_l, mask)
set_initial_cache! = make_set_initial_cache(model)
exp_tendency! = make_exp_tendency(model);
imp_tendency! = ClimaLand.make_imp_tendency(model);
update_jacobian! = ClimaLand.make_update_jacobian(model);

set_initial_cache!(p, Y, t0)
stepper = CTS.ARS111()
norm_condition = CTS.MaximumError(FT(1e-8))
conv_checker = CTS.ConvergenceChecker(; norm_condition = norm_condition)
ode_algo = CTS.IMEXAlgorithm(
    stepper,
    CTS.NewtonsMethod(
        max_iters = 2,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        convergence_checker = conv_checker,
    ),
)

# set up jacobian info
jac_kwargs =
    (; jac_prototype = RichardsTridiagonalW(Y), Wfact = update_jacobian!)

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
saveat = Array(t0:dt:tf)
sv = (;
    t = Array{Float64}(undef, length(saveat)),
    saveval = Array{NamedTuple}(undef, length(saveat)),
)
saving_cb = ClimaLand.NonInterpSavingCallback(sv, saveat)
updateat = Array(t0:dt:tf)
updatefunc = ClimaLand.make_update_drivers(atmos, nothing)
driver_cb = ClimaLand.DriverUpdateCallback(updateat, updatefunc)
cb = SciMLBase.CallbackSet(driver_cb, saving_cb)
@time sol = SciMLBase.solve(prob, ode_algo; dt = dt, saveat = dt, callback = cb)

longpts = range(-180.0, 180.0, 101)
latpts = range(-90.0, 90.0, 101)
hcoords = [
    ClimaCore.Geometry.LatLongPoint(lat, long) for long in longpts,
    lat in latpts
]
remapper = ClimaCore.Remapping.Remapper(surface_space, hcoords)

h∇_end = ClimaCore.Remapping.interpolate(
    remapper,
    oceans_to_zero.(sv.saveval[end].soil.h∇, mask),
)
fig = Figure(size = (600, 400))
ax = Axis(
    fig[1, 1],
    xlabel = "Longitude",
    ylabel = "Latitude",
    title = "Water table thickness",
)
clims = extrema(h∇_end)
CairoMakie.heatmap!(ax, longpts, latpts, h∇_end, colorrange = clims)
Colorbar(fig[:, end + 1], colorrange = clims)
outfile = joinpath(outdir, string("heatmap_h∇.png"))
CairoMakie.save(outfile, fig)

R_s_end = ClimaCore.Remapping.interpolate(
    remapper,
    oceans_to_zero.(sv.saveval[end].soil.R_s, mask),
)
fig = Figure(size = (600, 400))
ax = Axis(
    fig[1, 1],
    xlabel = "Longitude",
    ylabel = "Latitude",
    title = "Surface Runoff",
)
clims = extrema(R_s_end)
CairoMakie.heatmap!(ax, longpts, latpts, R_s_end, colorrange = clims)
Colorbar(fig[:, end + 1], colorrange = clims)
outfile = joinpath(outdir, string("heatmap_R_s.png"))
CairoMakie.save(outfile, fig)

R_ss_end = ClimaCore.Remapping.interpolate(
    remapper,
    oceans_to_zero.(sv.saveval[end].soil.R_ss, mask),
)

fig = Figure(size = (600, 400))
ax = Axis(
    fig[1, 1],
    xlabel = "Longitude",
    ylabel = "Latitude",
    title = "Subsurface Runoff",
)
clims = extrema(R_ss_end)
CairoMakie.heatmap!(ax, longpts, latpts, R_ss_end, colorrange = clims)
Colorbar(fig[:, end + 1], colorrange = clims)
outfile = joinpath(outdir, string("heatmap_R_ss.png"))
CairoMakie.save(outfile, fig)

θ_sfc_end = ClimaCore.Remapping.interpolate(
    remapper,
    ClimaLand.Soil.get_top_surface_field(
        oceans_to_zero.(sol.u[end].soil.ϑ_l, mask),
        surface_space,
    ),
)
fig = Figure(size = (1000, 400))
ax = Axis(fig[1, 1], xlabel = "Longitude", ylabel = "Latitude", title = "θ_sfc")
clims1 = extrema(θ_sfc_end)
CairoMakie.heatmap!(ax, longpts, latpts, θ_sfc_end, colorrange = clims1)
Colorbar(fig[1, 2], colorrange = clims1)

Δθ_sfc = ClimaCore.Remapping.interpolate(
    remapper,
    ClimaLand.Soil.get_top_surface_field(
        oceans_to_zero.(sol.u[end].soil.ϑ_l .- sol.u[1].soil.ϑ_l, mask),
        surface_space,
    ),
)
ax2 = Axis(fig[1, 3], xlabel = "Longitude", title = "θ_sfc Δ")
clims2 = extrema(Δθ_sfc)
CairoMakie.heatmap!(ax2, longpts, latpts, Δθ_sfc, colorrange = clims2)
Colorbar(fig[1, 4], colorrange = clims2)
outfile = joinpath(outdir, string("heatmap_θ_sfc.png"))
CairoMakie.save(outfile, fig)


int_ϑ_end = similar(sv.saveval[1].soil.h∇)
ClimaCore.Operators.column_integral_definite!(
    int_ϑ_end,
    oceans_to_zero.(sol.u[end].soil.ϑ_l, mask),
)
normalized_int_ϑ_end =
    ClimaCore.Remapping.interpolate(remapper, int_ϑ_end ./ 50.0)
fig = Figure(size = (1000, 400))
ax = Axis(
    fig[1, 1],
    xlabel = "Longitude",
    ylabel = "Latitude",
    title = "∫θ dz/Δz",
)
clims1 = extrema(normalized_int_ϑ_end)
CairoMakie.heatmap!(
    ax,
    longpts,
    latpts,
    normalized_int_ϑ_end,
    colorrange = clims1,
)
Colorbar(fig[1, 2], colorrange = clims1)

int_ϑ_1 = similar(sv.saveval[1].soil.h∇)
ClimaCore.Operators.column_integral_definite!(
    int_ϑ_1,
    oceans_to_zero.(sol.u[1].soil.ϑ_l, mask),
)
Δ_normalized_int_ϑ =
    ClimaCore.Remapping.interpolate(remapper, (int_ϑ_end .- int_ϑ_1) ./ 50.0)
ax2 = Axis(fig[1, 3], xlabel = "Longitude", title = "Δ∫θ dz /Δz")
clims2 = extrema(Δ_normalized_int_ϑ)
CairoMakie.heatmap!(
    ax2,
    longpts,
    latpts,
    Δ_normalized_int_ϑ,
    colorrange = clims2,
)
Colorbar(fig[1, 4], colorrange = clims2)
outfile = joinpath(outdir, string("heatmap_∫ϑdz.png"))
CairoMakie.save(outfile, fig)

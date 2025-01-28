# This shows how to run single column soil model, in standalone mode
# with spatially varying properties. We are mimicking the experiment
# carried out in Huang et. al.
# Can. J. Soil Sci. (2011) 91: 169183 doi:10.4141/CJSS09118,
# which measured the infiltration of layered soil in Fort McMurray,
# Alberta, Canada. We thank Mingbin Huang and S. Lee Barbour for
# correspondence and support, including sharing of data, with us.
# Note that all data used in this tutorial is available in their
# publication or included on this branch.

using CairoMakie
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import SciMLBase
import ClimaTimeSteppers as CTS
using ClimaCore
import ClimaParams as CP
using DelimitedFiles: readdlm

using ClimaLand
import ClimaLand.Simulations: LandSimulation, solve!
using ClimaLand.Domains: Column
using ClimaLand.Soil
import ClimaLand
FT = Float64;

# Define simulation times
t0 = Float64(0)
tf = Float64(60 * 60)
dt = Float64(30);

# Define the domain
zmax = FT(0)
zmin = FT(-1.1)
nelems = 50
Δ = FT((zmax - zmin) / nelems / 2)
soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems);

# Download the parameter data. This has been obtained from
# Table 1b of
# Infiltration and drainage processes in multi-layered coarse soils
# Mingbin Huang et. al.
# Can. J. Soil Sci. (2011) 91: 169183 doi:10.4141/CJSS09118
data_file = ClimaLand.Artifacts.huang_et_al2011_soil_van_genuchten_data();
parameter_data = readdlm(data_file, ',');
# Our model treats z as increasing in the upwards direction.
# Values below the surface are negative.
# Because of this, we convert the (positive-valued) depth
# of the data into a monotonically increasing z coordinate value.
# using a negative sign and the reverse function.
depth = reverse(-parameter_data[1, :] .* 0.01) # convert to m
ksat = reverse(parameter_data[6, :] .* 1 / 100.0 / 60.0) # convert cm/min to m/s
vgα = reverse(parameter_data[4, :] .* 100 * 2) # they report αᵈ; αʷ = 2αᵈ. This experiment is for infiltration (wetting).
vgn = reverse(parameter_data[5, :])
residual_frac = reverse(parameter_data[2, :])
porosity = reverse(parameter_data[3, :]);

# Create fields corresponding to the parameter
ν = SpaceVaryingInput(depth, porosity, soil_domain.space.subsurface)
K_sat = SpaceVaryingInput(depth, ksat, soil_domain.space.subsurface)
θ_r = SpaceVaryingInput(depth, residual_frac, soil_domain.space.subsurface);
# The specific storativity is not something we have data on, so we
# approximate it as being constant in depth, and create the parameter field directly:
S_s = ClimaCore.Fields.zeros(soil_domain.space.subsurface) .+ 1e-3;

# The retention model is a vanGenuchten model with α and n as a function of
# depth, read from the data:
hcm = SpaceVaryingInput(
    depth,
    (; α = vgα, n = vgn),
    soil_domain.space.subsurface,
    vanGenuchten{FT},
);

# The parameter struct:
params = ClimaLand.Soil.RichardsParameters(;
    ν = ν,
    hydrology_cm = hcm,
    K_sat = K_sat,
    S_s = S_s,
    θ_r = θ_r,
);

# From here on out, everything should look familiar if you've already gone
# through the other soil tutorials.
# Set Boundary conditions:
# At the top, we set moisture to equivalent of 5cm hydraulic head
# The bottom is free drainage
top_bc = ClimaLand.Soil.MoistureStateBC((p, t) -> 0.46705)
bottom_bc = ClimaLand.Soil.FreeDrainage()
boundary_fluxes = (; top = top_bc, bottom = bottom_bc)
soil = Soil.RichardsModel{FT}(;
    parameters = params,
    domain = soil_domain,
    boundary_conditions = boundary_fluxes,
    sources = (),
);

# Initial conditions read from data
data = readdlm("./docs/src/tutorials/standalone/Soil/sv_62_measurements.csv", ',')
data_z = Float64.(-1 .* data[2:end, 1]) ./ 100
data_ic = Float64.(data[2:end, 2])
data_8min = Float64.(data[2:end, 3])
data_16min = Float64.(data[2:end, 4])
data_24min = Float64.(data[2:end, 5])
data_32min = Float64.(data[2:end, 6])
data_40min = Float64.(data[2:end, 7])
data_60min = Float64.(data[2:end, 8])
function set_ic!(Y,p,t0,model; data_z, data_ic)
    Y.soil.ϑ_l .= SpaceVaryingInput(
        reverse(data_z),
        reverse(data_ic),
        axes(Y.soil.ϑ_l))
end
set_ic!(Y,p,t0,model) = set_ic!(Y,p,t0,model;data_z, data_ic)

# Timestepping:
stepper = CTS.ARS111()
err = 1e-8
convergence_cond = CTS.MaximumError(err)
conv_checker = CTS.ConvergenceChecker(norm_condition = convergence_cond)
ode_algo = CTS.IMEXAlgorithm(
    stepper,
    CTS.NewtonsMethod(
        max_iters = 10,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        convergence_checker = conv_checker,
    ),
)
saveat = [0.0, 8.0, 16.0, 24.0, 32.0, 40.0, 60.0] .* 60 # chosen to compare with data in plots in paper
simulation = LandSimulation(
    t0,
    tf,
    dt,
    soil;
    set_ic! = set_ic!,
    solver_kwargs = (; saveat),
    timestepper = ode_algo,
    user_callbacks = (),
    diagnostics = (),
);
sol = solve!(simulation);

z = parent(ClimaCore.Fields.coordinate_field(soil_domain.space.subsurface).z)[:]
ϑ_l = [parent(sol.u[k].soil.ϑ_l) for k in 1:length(sol.t)]

# Hydrus data

hydrus = readdlm("./docs/src/tutorials/standalone/Soil/sv_62_hydrus.csv", ',')
hydrus_z = Float64.(-1 .* hydrus[2:end, 1]) ./ 100
hydrus_ic = Float64.(hydrus[2:end, 2])
hydrus_8min = Float64.(hydrus[2:end, 3])
hydrus_16min = Float64.(hydrus[2:end, 4])
hydrus_24min = Float64.(hydrus[2:end, 5])
hydrus_32min = Float64.(hydrus[2:end, 6])
hydrus_40min = Float64.(hydrus[2:end, 7])
hydrus_60min = Float64.(hydrus[2:end, 8])


fig = CairoMakie.Figure(size = (600, 800), fontsize = 30)
ax = CairoMakie.Axis(
    fig[1, 1],
    xlabel = "Volumetric Water Content",
    ylabel = "Depth (m)",
    xgridvisible = false,
    ygridvisible = false,
    yticks = [-1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1],
)
lines!(ax, ϑ_l[1][:], z, label = "ClimaLand", color = "black", linewidth=3)
lines!(ax, ϑ_l[2][:], z, color = "blue", linewidth=3)
lines!(ax, ϑ_l[3][:], z, color = "green", linewidth=3)
lines!(ax, ϑ_l[4][:], z, color = "orange", linewidth=3)
lines!(ax, ϑ_l[5][:], z, color = "purple", linewidth=3)
lines!(ax, ϑ_l[6][:], z, color = "cyan", linewidth=3)
lines!(ax, ϑ_l[7][:], z, color = "brown", linewidth=3)

lines!(
    ax,
    hydrus_ic,
    hydrus_z,
    label = "Hydrus",
    color = "black",
    linestyle = :dash,
)
lines!(ax, hydrus_8min, hydrus_z, color = "blue", linestyle = :dash, linewidth=3)
lines!(ax, hydrus_16min, hydrus_z, color = "green", linestyle = :dash, linewidth=3)
lines!(ax, hydrus_24min, hydrus_z, color = "orange", linestyle = :dash, linewidth=3)
lines!(ax, hydrus_32min, hydrus_z, color = "purple", linestyle = :dash, linewidth=3)
lines!(ax, hydrus_40min, hydrus_z, color = "cyan", linestyle = :dash, linewidth=3)
lines!(ax, hydrus_60min, hydrus_z, color = "brown", linestyle = :dash, linewidth=3)

scatter!(
    ax,
    data_ic,
    data_z,
    label = "0 min",
    color = "black",
    marker = :xcross,
)
scatter!(
    ax,
    data_8min,
    data_z,
    color = "blue",
    label = "8 min",
    marker = :xcross,
)
scatter!(
    ax,
    data_16min,
    data_z,
    color = "green",
    label = "16 min",
    marker = :xcross,
)
scatter!(
    ax,
    data_24min,
    data_z,
    color = "orange",
    label = "24 min",
    marker = :xcross,
)
scatter!(
    ax,
    data_32min,
    data_z,
    color = "purple",
    label = "32 min",
    marker = :xcross,
)
scatter!(
    ax,
    data_40min,
    data_z,
    color = "cyan",
    label = "40 min",
    marker = :xcross,
)
scatter!(
    ax,
    data_60min,
    data_z,
    color = "brown",
    label = "60 min",
    marker = :xcross,
)
axislegend(ax, position = :rb, framevisible=false)
limits!(ax, 0, 0.7, -1.1, 0)

CairoMakie.save("./sv62_alpha_2_inf_updated_data_climaland_moisturebc.png", fig)
# ![](sv62_alpha_2_inf_updated_data_climaland_moisturebc.png)

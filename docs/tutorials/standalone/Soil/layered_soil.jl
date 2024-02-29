# This shows how to run single colum soil model, in standalone mode
# with spatially varying properties. We are mimicking the experiment
# carried out in Huang et. al.
# Can. J. Soil Sci. (2011) 91: 169183 doi:10.4141/CJSS09118,
# which measured the infiltration of layered soil in Fort McMurray,
# Alberta, Canada. We thank Mingbin Huang and S. Lee Barbour for
# correspondence and support, including sharing of data, with us
# Note that all data used in this tutorial is available in their
# publication.

using Plots
import SciMLBase
import ClimaTimeSteppers as CTS
using ClimaCore
import CLIMAParameters as CP
using ArtifactWrappers
using DelimitedFiles: readdlm

using ClimaLand
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
nelems = 75
Δ = FT((zmax - zmin) / nelems / 2)
soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems);

# Download the parameter data. This has been obtained from
# Table 1b of
# Infiltration and drainage processes in multi-layered coarse soils
# Mingbin Huang et. al.
# Can. J. Soil Sci. (2011) 91: 169183 doi:10.4141/CJSS09118
af = ArtifactFile(
    url = "https://caltech.box.com/shared/static/qvbt37xkwz8gveyi6tzzbs0e18trpcsq.csv",
    filename = "sv_62.csv",
)
dataset = ArtifactWrapper(@__DIR__, "sv62", ArtifactFile[af]);
dataset_path = get_data_folder(dataset);
data = joinpath(dataset_path, af.filename)
parameter_data = readdlm(data, ',');
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
# At the top, we use the observed value of Ksat at the top of the domain.
# Setting the flux to be -Ksat is approximating the top as saturated.
function top_flux_function(p, t)
    return -0.0001033
end
top_bc = ClimaLand.Soil.FluxBC(top_flux_function)
bottom_bc = ClimaLand.Soil.FreeDrainage()
boundary_fluxes = (; top = top_bc, bottom = bottom_bc)
soil = Soil.RichardsModel{FT}(;
    parameters = params,
    domain = soil_domain,
    boundary_conditions = boundary_fluxes,
    sources = (),
);

# Initial the state vectors, and set initial conditions
Y, p, cds = initialize(soil);

# Initial conditions
Y.soil.ϑ_l .= 0.0353; # read from Figure 4 of Huang et al. 

# We also set the initial conditions of the auxiliary state here:
set_initial_cache! = make_set_initial_cache(soil)
set_initial_cache!(p, Y, t0);

# Timestepping:
stepper = CTS.ARS111()
@assert FT in (Float32, Float64)
err = (FT == Float64) ? 1e-8 : 1e-4
convergence_cond = CTS.MaximumError(err)
conv_checker = CTS.ConvergenceChecker(norm_condition = convergence_cond)
ode_algo = CTS.IMEXAlgorithm(
    stepper,
    CTS.NewtonsMethod(
        max_iters = 50,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
        convergence_checker = conv_checker,
    ),
)
exp_tendency! = make_exp_tendency(soil)
imp_tendency! = make_imp_tendency(soil)
update_jacobian! = make_update_jacobian(soil)
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
saveat = [0.0, 8.0, 16.0, 24.0, 32.0, 40.0, 60.0] .* 60 # chosen to compare with data in plots in paper
sol = SciMLBase.solve(prob, ode_algo; dt = dt, saveat = saveat);

z = parent(ClimaCore.Fields.coordinate_field(soil_domain.space.subsurface).z)
ϑ_l = [parent(sol.u[k].soil.ϑ_l) for k in 1:length(sol.t)]
plot(ϑ_l[1], z, label = "initial", color = "grey", aspect_ratio = 0.8)
plot!(ϑ_l[2], z, label = "8min", color = "orange")
plot!(ϑ_l[3], z, label = "16min", color = "red")
plot!(ϑ_l[4], z, label = "24min", color = "teal")
plot!(ϑ_l[5], z, label = "32min", color = "blue")
plot!(ϑ_l[6], z, label = "40min", color = "purple")
plot!(ϑ_l[7], z, label = "60min", color = "green")
scatter!(porosity, depth, label = "Porosity")
plot!(legend = :bottomright)

plot!(xlim = [0, 0.7])

plot!(
    ylim = [-1.1, 0],
    yticks = [-1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1],
)

plot!(ylabel = "Depth (m)")

plot!(xlabel = "Volumeteric Water Content")

savefig("./sv62_alpha_2_inf_updated_data_climaland.png")
# ![](sv62_alpha_2_inf_updated_data_climaland.png)

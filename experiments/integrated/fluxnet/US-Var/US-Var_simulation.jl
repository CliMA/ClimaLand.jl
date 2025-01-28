"""This file contains simulation variables for running Clima Land on the US-Var
fluxtower site. This includes both the domain variables and timestepping
variables for running the simulation."""

# DOMAIN SETUP:

# Column dimensions
# Stem and leaf compartments and their heights:
n_stem = Int64(0)
n_leaf = Int64(1)
h_leaf = FT(0.5) # m, Xu and Baldocchi, 2003
h_stem = FT(0) # m

# TIME STEPPING:

# Starting time:
t0 = Float64(0)

# Time step size:
dt = Float64(900)
# For soil column
nelements = 14
zmin = FT(-0.5) #m, Xu and Baldocchi, 2003
zmax = FT(0)
h_canopy = h_stem + h_leaf
compartment_midpoints =
    n_stem > 0 ? [h_stem / 2, h_stem + h_leaf / 2] : [h_leaf / 2]
compartment_surfaces = n_stem > 0 ? [zmax, h_stem, h_canopy] : [zmax, h_leaf]

land_domain = Column(;
    zlim = (zmin, zmax),
    nelements = nelements,
)
canopy_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)

N_spinup_days = 30
N_days = N_spinup_days + 365

tf = Float64(t0 + 3600 * 24 * N_days)
t_spinup = Float64(t0 + N_spinup_days * 3600 * 24)

# Set up timestepper
timestepper = CTS.ARS111();
ode_algo = CTS.IMEXAlgorithm(
    timestepper,
    CTS.NewtonsMethod(
        max_iters = 6,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
    ),
);

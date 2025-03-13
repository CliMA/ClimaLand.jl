using ClimaComms
ClimaComms.@import_required_backends
import SciMLBase
import ClimaTimeSteppers as CTS
using ClimaCore.MatrixFields
import ClimaCore.MatrixFields: @name
using Dates
using Test
import ClimaParams as CP
using ClimaLand
import ClimaLand.Parameters as LP

context = ClimaComms.context()
nelements = (101, 15)
start_date = DateTime(2008)
Δt = 450.0
t0 = 0.0
tf = t0 + Δt

FT = Float64
earth_param_set = LP.LandParameters(FT)

f_over = FT(3.28) # 1/m
R_sb = FT(1.484e-4 / 1000) # m/s
other_soil_params = (; f_over, R_sb)

α_snow = FT(0.67)
other_snow_params = (; α_snow,)

# Energy Balance model
ac_canopy = FT(2.5e3)
# Plant Hydraulics and general plant parameters
SAI = FT(0.0) # m2/m2
f_root_to_shoot = FT(3.5)
RAI = FT(1.0)
K_sat_plant = FT(5e-9) # m/s # seems much too small?
ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value is -4 MPa
Weibull_param = FT(4) # unitless, Holtzman's original c param value
a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
plant_ν = FT(1.44e-4)
plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
n_stem = 0
n_leaf = 1
h_stem = FT(0.0)
h_leaf = FT(1.0)
zmax = FT(0.0)
other_canopy_params = (;
    ac_canopy,
    SAI,
    f_root_to_shoot,
    RAI,
    K_sat_plant,
    a,
    ψ63,
    Weibull_param,
    plant_ν,
    plant_S_s,
    n_stem,
    n_leaf,
    h_stem,
    h_leaf,
    zmax,
)

land, Y, p, cds = ClimaLand.land_model_setup(
    FT;
    earth_param_set,
    context,
    nelements,
    start_date,
    t0,
    Δt,
    other_soil_params,
    other_canopy_params,
    other_snow_params,
    era5_lowres = true,
)

surface_space = axes(Y.snow.U)
binary_mask = .~parent(surface_space.grid.mask.is_active)[:]


# First, set the cache and make sure there are no NaNs, or values set over the ocean
set_initial_cache! = make_set_initial_cache(land)
set_initial_cache!(p, Y, t0)

for var in propertynames(p.soil)
    field_values = parent(getproperty(p.soil, var))
    @test sum(isnan, field_values) == 0
    if length(size(field_values)) == 5 # 3d var
        @test extrema(field_values[:, 1, 1, :, binary_mask]) == (0.0, 0.0)
    else
        @test extrema(field_values[1, 1, :, binary_mask]) == (0.0, 0.0)
    end
end
for var in propertynames(p.soilco2)
    field_values = parent(getproperty(p.soilco2, var))
    @test sum(isnan, field_values) == 0
    if length(size(field_values)) == 5 # 3d var
        @test extrema(field_values[:, 1, 1, :, binary_mask]) == (0.0, 0.0)
    else
        @test extrema(field_values[1, 1, :, binary_mask]) == (0.0, 0.0)
    end
end
for var in propertynames(p.snow)
    field_values = parent(getproperty(p.snow, var))
    @test sum(isnan, field_values) == 0
    @test extrema(field_values[1, 1, 1, binary_mask]) == (0.0, 0.0)
end
for var in propertynames(p.drivers)
    field_values = parent(getproperty(p.drivers, var))
    #@test extrema(field_values[1,1,1,binary_mask]) == (0.0,0.0) # FAILS
    @test sum(isnan, field_values) == 0
end
field_pn_p = [
    pn for pn in propertynames(p) if pn != :soil &&
    pn != :canopy &&
    pn != :snow &&
    pn != :soilco2 &&
    pn != :drivers &&
    ~occursin("dss", String(pn))
]
for var in field_pn_p
    field_values = parent(getproperty(p, var))
    @test sum(isnan, field_values) == 0
    if length(size(field_values)) == 5 # 3d var
        @test extrema(field_values[:, 1, 1, :, binary_mask]) == (0.0, 0.0)
    else
        @test extrema(field_values[1, 1, :, binary_mask]) == (0.0, 0.0)
    end
end
for key in keys(p.canopy)
    x = getproperty(p.canopy, key)
    for var in propertynames(x)
        field_values = parent(getproperty(x, var))
        @test extrema(field_values[1, 1, :, binary_mask]) == (0.0, 0.0)
        @test sum(isnan, field_values) == 0
    end
end


# now check what tendency updates do:
dY = similar(Y)
@. dY = 0

imp_tendency! = make_imp_tendency(land)
imp_tendency!(dY, Y, p, t0)
# Test that the masked parts of dY did not update and are still zero
@test extrema(parent(dY.soil.ϑ_l)[:, 1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(dY.soil.θ_i)[:, 1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(dY.soil.ρe_int)[:, 1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(dY.snow.U)[1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(dY.snow.S)[1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(dY.snow.S_l)[1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(dY.canopy.energy.T)[1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(
    Array(parent(dY.canopy.hydraulics.ϑ_l.:1))[1, 1, 1, Array(binary_mask)],
) == (0.0, 0.0)

@. dY = 0
exp_tendency! = make_exp_tendency(land)
exp_tendency!(dY, Y, p, t0)
# Test that the masked parts of dY did not update and are still zero
@test extrema(parent(dY.soil.ϑ_l)[:, 1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(dY.soil.θ_i)[:, 1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(dY.soil.ρe_int)[:, 1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(dY.snow.U)[1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(dY.snow.S)[1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(dY.snow.S_l)[1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(dY.canopy.energy.T)[1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(
    Array(parent(dY.canopy.hydraulics.ϑ_l.:1))[1, 1, 1, Array(binary_mask)],
) == (0.0, 0.0)

# Check jacobian

jacobian! = ClimaLand.make_jacobian(land);
jac_prototype = ClimaLand.FieldMatrixWithSolver(Y);
# Check the cache fields
@test axes(
    jac_prototype.solver.cache.cache₁.entries[1],
).grid.horizontal_grid.mask == surface_space.grid.mask
for i in 1:3
    @test axes(
        jac_prototype.solver.cache.b₂′.entries[i],
    ).grid.horizontal_grid.mask == surface_space.grid.mask
    @test axes(
        jac_prototype.solver.cache.cache₂.entries[i],
    ).grid.horizontal_grid.mask == surface_space.grid.mask
end
for i in 4:8
    @test axes(jac_prototype.solver.cache.b₂′.entries[i]).grid.mask ==
          surface_space.grid.mask
    @test axes(jac_prototype.solver.cache.cache₂.entries[i]).grid.mask ==
          surface_space.grid.mask

end

jacobian!(jac_prototype, Y, p, Δt, t0);
(; matrix) = jac_prototype;
∂ϑres∂ϑ = matrix[@name(soil.ϑ_l), @name(soil.ϑ_l)];
@test extrema(
    Array(parent(∂ϑres∂ϑ.entries.:1))[:, 1, 1, 1, Array(binary_mask)],
) == (0.0, 0.0)
@test extrema(
    Array(parent(∂ϑres∂ϑ.entries.:2))[:, 1, 1, 1, Array(binary_mask)],
) == (0.0, 0.0)
@test extrema(
    Array(parent(∂ϑres∂ϑ.entries.:3))[:, 1, 1, 1, Array(binary_mask)],
) == (0.0, 0.0)

∂ρeres∂ρe = matrix[@name(soil.ρe_int), @name(soil.ρe_int)];
@test extrema(
    Array(parent(∂ρeres∂ρe.entries.:1))[:, 1, 1, 1, Array(binary_mask)],
) == (0.0, 0.0)
@test extrema(
    Array(parent(∂ρeres∂ρe.entries.:2))[:, 1, 1, 1, Array(binary_mask)],
) == (0.0, 0.0)
@test extrema(
    Array(parent(∂ρeres∂ρe.entries.:3))[:, 1, 1, 1, Array(binary_mask)],
) == (0.0, 0.0)
∂ρeres∂ϑ = matrix[@name(soil.ρe_int), @name(soil.ϑ_l)];
@test extrema(
    Array(parent(∂ρeres∂ϑ.entries.:1))[:, 1, 1, 1, Array(binary_mask)],
) == (0.0, 0.0)
@test extrema(
    Array(parent(∂ρeres∂ϑ.entries.:2))[:, 1, 1, 1, Array(binary_mask)],
) == (0.0, 0.0)
@test extrema(
    Array(parent(∂ρeres∂ϑ.entries.:3))[:, 1, 1, 1, Array(binary_mask)],
) == (0.0, 0.0)

∂Tres∂T = matrix[@name(canopy.energy.T), @name(canopy.energy.T)];
@test extrema(Array(parent(∂Tres∂T.entries.:1))[1, 1, 1, Array(binary_mask)]) ==
      (0.0, 0.0)



b = similar(Y) .* 1;
x = deepcopy(Y);
@test axes(x.soil.ϑ_l).grid.horizontal_grid.mask == surface_space.grid.mask;


MatrixFields.field_matrix_solve!(
    jac_prototype.solver,
    x,
    jac_prototype.matrix,
    b,
);

@test extrema(parent(x.soil.ϑ_l)[:, 1, 1, 1, binary_mask]) == (0.0, 0.0) # FAILS
@test extrema(parent(x.soil.θ_i)[:, 1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(x.soil.ρe_int)[:, 1, 1, 1, binary_mask]) == (0.0, 0.0) # FAILS
@test extrema(parent(x.snow.U)[1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(x.snow.S)[1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(x.snow.S_l)[1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(x.canopy.energy.T)[1, 1, 1, binary_mask]) == (0.0, 0.0) #FAILS
@test extrema(
    Array(parent(x.canopy.hydraulics.ϑ_l.:1))[1, 1, 1, Array(binary_mask)],
) == (0.0, 0.0)

# Take a step
jac_kwargs = (; jac_prototype = jac_prototype, Wfact = jacobian!)
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
# Define timestepper and ODE algorithm
stepper = CTS.ARS111()
ode_algo = CTS.IMEXAlgorithm(
    stepper,
    CTS.NewtonsMethod(
        max_iters = 3,
        update_j = CTS.UpdateEvery(CTS.NewNewtonIteration),
    ),
)

sol = SciMLBase.solve(
    prob,
    ode_algo;
    dt = Δt,
    adaptive = false,
    saveat = [t0, tf],
);
u = sol.u[end];
@test extrema(parent(u.soil.ϑ_l)[:, 1, 1, 1, binary_mask]) == (0.0, 0.0) # FAILS
@test extrema(parent(u.soil.θ_i)[:, 1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(u.soil.ρe_int)[:, 1, 1, 1, binary_mask]) == (0.0, 0.0) # FAILS
@test extrema(parent(u.snow.U)[1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(u.snow.S)[1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(u.snow.S_l)[1, 1, 1, binary_mask]) == (0.0, 0.0)
@test extrema(parent(u.canopy.energy.T)[1, 1, 1, binary_mask]) == (0.0, 0.0) # FAILS
@test extrema(
    Array(parent(u.canopy.hydraulics.ϑ_l.:1))[1, 1, 1, Array(binary_mask)],
) == (0.0, 0.0)

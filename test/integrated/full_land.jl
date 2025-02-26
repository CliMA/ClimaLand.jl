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
using ClimaCore
include("full_land_setup.jl")
using ClimaUtilities.ClimaArtifacts

import Interpolations
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.Regridders: InterpolationsRegridder
using ClimaCore: Spaces

function test_p(p, binary_mask)
    for var in propertynames(p.drivers)
        @show var
        field_values = parent(getproperty(p.drivers, var))
        if var != :soc
            @test extrema(field_values[1, 1, 1, binary_mask]) == (0.0, 0.0)
        else
            @test extrema(field_values[:, 1, 1, 1, binary_mask]) == (0.0, 0.0)
        end
        @test sum(isnan, field_values) == 0
    end

    for var in propertynames(p.soil)
        @show var
        field_values = parent(getproperty(p.soil, var))
        @test sum(isnan, field_values) == 0
        if length(size(field_values)) == 5 # 3d var
            @test extrema(field_values[:, 1, 1, :, binary_mask]) == (0.0, 0.0)
        else
            @test extrema(field_values[1, 1, :, binary_mask]) == (0.0, 0.0)
        end
    end
    for var in propertynames(p.soilco2)
        @show var
        field_values = parent(getproperty(p.soilco2, var))
        @test sum(isnan, field_values) == 0
        if length(size(field_values)) == 5 # 3d var
            @test extrema(field_values[:, 1, 1, :, binary_mask]) == (0.0, 0.0)
        else
            @test extrema(field_values[1, 1, :, binary_mask]) == (0.0, 0.0)
        end
    end
    for var in propertynames(p.snow)
        @show var
        field_values = parent(getproperty(p.snow, var))
        @test sum(isnan, field_values) == 0
        @test extrema(field_values[1, 1, 1, binary_mask]) == (0.0, 0.0)
    end

    field_pn_p = [
        pn for pn in propertynames(p) if pn != :soil &&
        pn != :canopy &&
        pn != :snow &&
        pn != :soilco2 &&
        pn != :drivers &&
        ~occursin("dss", String(pn))
    ]
    for key in keys(p.canopy)
        @show key
        x = getproperty(p.canopy, key)
        for var in propertynames(x)
            field_values = parent(getproperty(x, var))
            @test extrema(field_values[1, 1, :, binary_mask]) == (0.0, 0.0)
            @test sum(isnan, field_values) == 0
        end
    end
    for var in field_pn_p
        @show var
        field_values = parent(getproperty(p, var))
        @test sum(isnan, field_values) == 0
        if length(size(field_values)) == 5 # 3d var
            @test extrema(field_values[:, 1, 1, :, binary_mask]) == (0.0, 0.0)
        else
            @test extrema(field_values[1, 1, :, binary_mask]) == (0.0, 0.0)
        end
    end

end


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
scalar_soil_params = (; f_over, R_sb)

α_snow = FT(0.67)
scalar_snow_params = (; α_snow, Δt)

# Energy Balance model
ac_canopy = FT(2.5e3)
SAI = FT(0.0) # m2/m2
f_root_to_shoot = FT(3.5)
RAI = FT(1.0)
K_sat_plant = FT(5e-9) # m/s # seems much too small?
ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value is -4 MPa
Weibull_param = FT(4) # unitless, Holtzman's original c param value
a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
plant_ν = FT(1.44e-4)
plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
h_leaf = FT(1.0)
n_stem = 0
n_leaf = 1
h_stem = FT(0.0)
h_leaf = FT(1.0)
zmax = FT(0.0)

scalar_canopy_params = (;
    ac_canopy,
    K_sat_plant,
    a,
    ψ63,
    Weibull_param,
    plant_ν,
    plant_S_s,
    h_leaf,
)

domain = ClimaLand.global_domain(FT; nelements = nelements)
surface_space = domain.space.surface
start_date = DateTime(2008)
# Forcing data
era5_ncdata_path = joinpath(
    ClimaLand.Artifacts.era5_land_forcing_data2008_folder_path(;
        context,
        lowres = true,
    ),
    "era5_2008_1.0x1.0_lowres.nc",
)
forcing = ClimaLand.prescribed_forcing_era5(
    joinpath(era5_ncdata_path),
    surface_space,
    start_date,
    earth_param_set,
    FT,
)
LAI = ClimaLand.prescribed_lai_modis(
    joinpath(
        ClimaLand.Artifacts.modis_lai_forcing_data_path(; context),
        "Yuan_et_al_2008_1x1.nc",
    ),
    domain.space.surface,
    start_date,
)

land = global_land_model(
    FT,
    scalar_soil_params,
    scalar_canopy_params,
    scalar_snow_params,
    earth_param_set;
    context = context,
    domain = domain,
    forcing = forcing,
    LAI = LAI,
)

Y, p, cds = initialize(land)

@testset "Total energy and water" begin
    # Soil IC
    ϑ_l0 = land.soil.parameters.ν ./ 2
    θ_i0 = land.soil.parameters.ν ./ 5
    T = FT(270.0)
    ρc_s = @. Soil.volumetric_heat_capacity(
        ϑ_l0,
        θ_i0,
        land.soil.parameters.ρc_ds,
        earth_param_set,
    )
    ρe_int0 = @. Soil.volumetric_internal_energy(θ_i0, ρc_s, T, earth_param_set)
    Y.soil.ϑ_l .= ϑ_l0
    Y.soil.θ_i .= θ_i0

    Y.soil.ρe_int = ρe_int0
    
    # Canopy IC
    ϑ0 = land.canopy.hydraulics.parameters.ν / 2
    CTemp0 = FT(290.5)
    
    Y.canopy.hydraulics.ϑ_l.:1 .= ϑ0
    
    Y.canopy.energy.T .= CTemp0
    
    # Snow IC
    S0 = FT(0.5)
    S_l0 = FT(0.3)
    STemp0 = FT(270)
    U0 = Snow.energy_from_T_and_swe(S0, STemp0, land.snow.parameters)
    
    Y.snow.S .= S0
    Y.snow.S_l .= S_l0
    Y.snow.U .= U0
    
    t0 = 0.0
    set_initial_cache! = make_set_initial_cache(land)
    set_initial_cache!(p, Y, t0)
    
    # Check total
    area_index = p.canopy.hydraulics.area_index.leaf
    h_canopy = land.canopy.hydraulics.compartment_surfaces[end]
    ρ_ice = LP.ρ_cloud_ice(earth_param_set)
    ρ_liq = LP.ρ_cloud_liq(earth_param_set)
    int_cache = ClimaCore.Fields.zeros(domain.space.surface)
    ClimaCore.Operators.column_integral_definite!(
        int_cache,
        @. (ϑ_l0 + θ_i0 * ρ_ice / ρ_liq)
    )
    soil_exp = int_cache
    canopy_exp = @. area_index * h_canopy .* ϑ0
    snow_exp = S0
    total_water = ClimaCore.Fields.zeros(domain.space.surface)
    cache = ClimaCore.Fields.zeros(domain.space.surface)
    ClimaLand.total_liq_water_vol_per_area!(total_water, land, Y, p, t0, cache)
    @test all(parent(total_water) .≈ parent(snow_exp .+ canopy_exp .+ soil_exp))
    
    int_cache .*= 0
    ClimaCore.Operators.column_integral_definite!(int_cache, ρe_int0)
    soil_exp = int_cache
    canopy_exp = @. area_index * land.canopy.energy.parameters.ac_canopy * CTemp0
    snow_exp = U0
    total_energy = ClimaCore.Fields.zeros(domain.space.surface)
    ClimaLand.total_energy_per_area!(total_energy, land, Y, p, t0, cache)
    @test all(parent(total_energy) .≈ parent(snow_exp .+ canopy_exp .+ soil_exp))
end


@testset "Mask of full land" begin
    Y, p, cds = initialize(land)
    # Y Here we change them to be -1 for tracking errors below
    Y .= 0
    surface_space = axes(Y.snow.U)
    subsurface_space = axes(Y.soil.ϑ_l)
    binary_mask = .~parent(surface_space.grid.mask.is_active)[:]
    test_p(p, binary_mask)
    #=
    ic_path = ClimaLand.Artifacts.soil_ic_2008_50m_path(; context = context)
    ClimaLand.set_soil_initial_conditions!(Y, land.soil.parameters.ν, land.soil.parameters.θ_r, subsurface_space, ic_path)
    evaluate!(p.snow.T, land.snow.boundary_conditions.atmos.T, t0)
    ClimaLand.set_snow_initial_conditions!(
        Y,
    p,
    surface_space,
        ic_path,
    land.snow.parameters,
    )
    
    Y.soilco2.C .= FT(0.000412) # set to atmospheric co2, mol co2 per mol air
    Y.canopy.hydraulics.ϑ_l.:1 .= plant_ν
    evaluate!(Y.canopy.energy.T, land.snow.boundary_conditions.atmos.T, t0)
    
    
    # First, set the cache and make sure there are no NaNs, or values set over the ocean
    set_initial_cache! = make_set_initial_cache(land)
    set_initial_cache!(p, Y, t0)
    test_p(p, binary_mask)
    
    
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
=#
end

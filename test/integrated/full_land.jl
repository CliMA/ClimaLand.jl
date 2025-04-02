using ClimaComms
ClimaComms.@import_required_backends
import SciMLBase
import ClimaTimeSteppers as CTS
using ClimaCore.MatrixFields
import ClimaCore.MatrixFields: @name
using ClimaUtilities.ClimaArtifacts
using Dates
using Test
import ClimaParams as CP
using ClimaLand
import ClimaLand.Parameters as LP
using ClimaCore
include("full_land_setup.jl");

#=
import ClimaCore: Spaces, Fields
import ClimaCore.DataLayouts: VIJFH, AbstractData
Base.map(fn, data::VIJFH...) =
	Base.map!(fn, similar(first(data)), data...)
function Base.map!(fn, dest::VIJFH, data::VIJFH...)
    (_, _, _, _, Nh) = size(first(data))
    @inbounds for h in 1:Nh, j in 1:Nij, i in 1:Nij, v in 1:Nv
        idx = CartesianIndex(i, j, 1, v, h)
        dest[idx] = fn(idx)
    end
    return dest
end
check_ocean(space::Spaces.AbstractSpace, field::Fields.Field) =
	check_ocean(space, Fields.field_values(field))
function check_ocean(space::Spaces.AbstractSpace, data::AbstractData)
	(; is_active) = Spaces.get_mask(space)
	map(is_active, data) do idx
		if !is_active[idx]
			@assert data[idx] == 0
		end
		data[idx] # need to return similar element
	end
end
=#
function test_p(p, binary_mask; val = 0.0)
    properties = [p.drivers, p.soil, p.soilco2, p.snow, p.canopy.energy, p.canopy.hydraulics, p.canopy.radiative_transfer, p.canopy.photosynthesis, p.canopy.sif,p.canopy.turbulent_fluxes, p.canopy.autotrophic_respiration, p.canopy.conductance]
    for property in properties
        for var in propertynames(property)
            field_values = parent(getproperty(property, var))
            if length(size(field_values)) == 5 # 3d var
                extrema(field_values[:, 1, 1, :, binary_mask]) == (val, val) ? nothing : @warn "$var not zero"
            else
                extrema(field_values[ 1, 1, :, binary_mask]) == (val, val) ? nothing : @warn "$var not zero"
            end
            @test sum(isnan, field_values) == 0
        end
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
        if length(size(field_values)) == 5 # 3d var
            extrema(field_values[:, 1, 1, :, binary_mask]) == (val, val) ? nothing : @warn "$var not zero"
        else
            extrema(field_values[ 1, 1, :, binary_mask]) == (val, val) ? nothing : @warn "$var not zero"
        end
        @test sum(isnan, field_values) == 0
    end
end

function test_field_vector(dY, binary_mask; val = 0.0)
     extrema(parent(dY.soil.ϑ_l)[:, 1, 1, 1, binary_mask]) == (val, val) ? nothing :  @warn "ϑ_l not zero"
     extrema(parent(dY.soil.θ_i)[:, 1, 1, 1, binary_mask]) == (val, val) ? nothing :  @warn "θ_i not zero"
     extrema(parent(dY.soil.ρe_int)[:, 1, 1, 1, binary_mask]) == (val, val) ? nothing :  @warn "ρe_int not zero"
     extrema(parent(dY.snow.U)[1, 1, 1, binary_mask]) == (val, val) ? nothing :  @warn "U not zero"
     extrema(parent(dY.snow.S)[1, 1, 1, binary_mask]) == (val, val) ? nothing :  @warn "S not zero"
     extrema(parent(dY.snow.S_l)[1, 1, 1, binary_mask]) == (val, val) ? nothing :  @warn "S_l not zero"
     extrema(parent(dY.canopy.energy.T)[1, 1, 1, binary_mask]) == (val, val) ? nothing :  @warn "Canopy T not zero"
     extrema(
        Array(parent(dY.canopy.hydraulics.ϑ_l.:1))[1, 1, 1, Array(binary_mask)],
    ) == (val, val) ? nothing :  @warn "Canopy ϑ_l not zero"
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
K_sat_plant = FT(5e-9) # m/s # seems much too small?
ψ63 = FT(-4 / 0.0098) # / MPa to m, Holtzman's original parameter value is -4 MPa
Weibull_param = FT(4) # unitless, Holtzman's original c param value
a = FT(0.05 * 0.0098) # Holtzman's original parameter for the bulk modulus of elasticity
plant_ν = FT(1.44e-4)
plant_S_s = FT(1e-2 * 0.0098) # m3/m3/MPa to m3/m3/m
h_leaf = FT(1.0)

scalar_canopy_params = (;
    ac_canopy,
    K_sat_plant,
    a,
    ψ63,
    Weibull_param,
    plant_ν,
    plant_S_s,
    h_leaf,
);

domain = ClimaLand.global_domain(FT; nelements = nelements);
surface_space = domain.space.surface;
start_date = DateTime(2008);
# Forcing data
era5_ncdata_path = joinpath(
    ClimaLand.Artifacts.era5_land_forcing_data2008_folder_path(;
        context,
        lowres = true,
    ),
    "era5_2008_1.0x1.0_lowres.nc",
);
forcing = ClimaLand.prescribed_forcing_era5(
    joinpath(era5_ncdata_path),
    surface_space,
    start_date,
    earth_param_set,
    FT,
);
LAI = ClimaLand.prescribed_lai_modis(
    joinpath(
        ClimaLand.Artifacts.modis_lai_forcing_data_path(; context),
        "Yuan_et_al_2008_1x1.nc",
    ),
    domain.space.surface,
    start_date,
);

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
);
#=
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

    oceans = .~parent(domain.space.surface.grid.mask.is_active)[:]
    continents = parent(domain.space.surface.grid.mask.is_active)[:]
    expected= snow_exp .+ canopy_exp .+ soil_exp # if we take parent first, this fails?
    @test all(parent(total_water)[1, 1, 1, continents] .≈ parent(expected)[1,1,1,continents])
    @test all(parent(total_water)[1, 1, 1, oceans] .≈ parent(expected)[1,1,1,oceans])
    
    int_cache .*= 0
    ClimaCore.Operators.column_integral_definite!(int_cache, ρe_int0)
    soil_exp = int_cache
    canopy_exp = @. area_index * land.canopy.energy.parameters.ac_canopy * CTemp0
    snow_exp = U0
    total_energy = ClimaCore.Fields.zeros(domain.space.surface)
    ClimaLand.total_energy_per_area!(total_energy, land, Y, p, t0, cache)
    expected = snow_exp .+ canopy_exp .+ soil_exp
    @test all(parent(total_energy)[1, 1, 1, continents] .≈ parent(expected)[1,1,1,continents])
    @test all(parent(total_energy)[1, 1, 1, oceans] .≈ parent(expected)[1,1,1,oceans])
end
=#

#@testset "Mask of full land" begin
    Y, p, cds = initialize(land);
    Y .= 0;
    surface_space = axes(Y.snow.U)
    subsurface_space = axes(Y.soil.ϑ_l)
    binary_mask = .~parent(surface_space.grid.mask.is_active)[:]
    # Test that the cache is zero over the ocean
    @info("testing initial cache")
    test_p(p, binary_mask)

    # Set initial conditions
    ic_path = ClimaLand.Artifacts.soil_ic_2008_50m_path(; context = context)
    T_bounds = (273.0, 290.0)
    ClimaLand.set_soil_initial_conditions!(Y, land.soil.parameters.ν, land.soil.parameters.θ_r, subsurface_space, ic_path, land.soil, T_bounds)
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
    
    
    # Now, set the cache with physical values and make sure there are no NaNs, or values set over the ocean
    set_initial_cache! = make_set_initial_cache(land)
    set_initial_cache!(p, Y, t0)
    @info("testing set cache")
    test_p(p, binary_mask)

    # Check tendency functions do not update the state over the ocean
    dY = similar(Y)    
    # Implicit tendency
    @. dY = 0
    imp_tendency! = make_imp_tendency(land)
    imp_tendency!(dY, Y, p, t0)
    @info("testing implicit tendency")
    test_field_vector(dY, binary_mask)

    # Explicit tendency
    @. dY = 0
    exp_tendency! = make_exp_tendency(land)
    exp_tendency!(dY, Y, p, t0)
    @info("testing explicit tendency")
    test_field_vector(dY, binary_mask)
    

    # Jacobian checks
    jacobian! = ClimaLand.make_jacobian(land);
    jac_prototype = ClimaLand.FieldMatrixWithSolver(Y);
    # Check that the jacobian cache fields inherit the mask
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

    # Check that the jacobian update respects the mask
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
    
    
    # Now carry out a solve of Jx = b, with x = Y, and b = 1
    b = similar(Y);
    fill!(parent(b), 1)
    x = deepcopy(Y);
    fill!(parent(x), 1)
    @test axes(x.soil.ϑ_l).grid.horizontal_grid.mask == surface_space.grid.mask;
    @test axes(b.soil.ϑ_l).grid.horizontal_grid.mask == surface_space.grid.mask;
    test_field_vector(x, binary_mask; val = 1.0)
    test_field_vector(b, binary_mask; val = 1.0)

    
    MatrixFields.field_matrix_solve!(jac_prototype.solver, x, jac_prototype.matrix, b);
    test_field_vector(x, binary_mask; val = 1.0)
    #=
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
    test_field_vector(u, binary_mask)
=#
#end

=#
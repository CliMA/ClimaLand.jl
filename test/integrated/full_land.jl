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
using ClimaLand.Domains
using ClimaLand.Soil
using ClimaLand.Canopy
using ClimaLand.Snow
import ClimaLand.Parameters as LP
using ClimaCore

for FT in (Float32, Float64)
    @testset "Default LandModel constructor, FT=$FT" begin
        toml_dict = LP.create_toml_dict(FT)
        domain = Domains.global_domain(FT)
        atmos, radiation = ClimaLand.prescribed_analytic_forcing(FT; toml_dict)
        forcing = (; atmos, radiation)
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2)

        # Soil model
        soil = Soil.EnergyHydrology{FT}(
            domain,
            forcing,
            toml_dict;
            prognostic_land_components,
            additional_sources = (ClimaLand.RootExtraction{FT}(),),
        )

        # SoilCO2 model
        co2_prognostic_soil =
            Soil.Biogeochemistry.PrognosticMet(soil.parameters)
        soil_organic_carbon =
            PrescribedSoilOrganicCarbon{FT}(TimeVaryingInput((t) -> 5))
        soilco2_drivers = Soil.Biogeochemistry.SoilDrivers(
            co2_prognostic_soil,
            soil_organic_carbon,
            atmos,
        )
        soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(
            domain,
            soilco2_drivers,
            toml_dict,
        )

        # Canopy model
        surface_domain = Domains.obtain_surface_domain(domain)
        LAI = TimeVaryingInput((t) -> FT(1.0))
        ground = ClimaLand.PrognosticGroundConditions{FT}()
        canopy_forcing = (; atmos, radiation, ground)
        canopy = Canopy.CanopyModel{FT}(
            surface_domain,
            canopy_forcing,
            LAI,
            toml_dict;
            prognostic_land_components,
        )

        # Snow model
        dt = FT(180)
        snow = SnowModel(
            FT,
            surface_domain,
            forcing,
            toml_dict,
            dt;
            prognostic_land_components,
        )

        model = LandModel{FT}(canopy, snow, soil, soilco2)

        # The constructor has many asserts that check the model
        # components, so we don't need to check them again here.
        @test model.soil == soil
        @test model.soilco2 == soilco2
        @test model.canopy == canopy
        @test model.snow == snow
    end
end

"""
    check_ocean_values_p(p, binary_mask; val = 0.0)

This function tests that every field stored in `p` has all of
its values (where binary_mask == 1) equal to  `val`. Note that
this is meant to be used with the full land model (canopy,
snow, soil, soilco2).

Useful for checking if land model functions are updating the
values over the ocean.
"""
function check_ocean_values_p(p, binary_mask; val = 0.0)
    properties = [
        p.drivers,
        p.soil,
        p.soilco2,
        p.snow,
        p.canopy.energy,
        p.canopy.hydraulics,
        p.canopy.radiative_transfer,
        p.canopy.photosynthesis,
        p.canopy.sif,
        p.canopy.turbulent_fluxes,
        p.canopy.autotrophic_respiration,
        p.canopy.conductance,
    ]
    for property in properties
        for var in propertynames(property)
            field_values = Array(parent(getproperty(property, var)))
            if length(size(field_values)) == 5 # 3d var
                @test extrema(field_values[:, 1, 1, :, Array(binary_mask)]) ==
                      (val, val)
            else
                @test extrema(field_values[1, 1, :, Array(binary_mask)]) ==
                      (val, val)
            end
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
        field_values = Array(parent(getproperty(p, var)))
        if length(size(field_values)) == 5 # 3d var
            @test extrema(field_values[:, 1, 1, :, Array(binary_mask)]) ==
                  (val, val)
        else
            @test extrema(field_values[1, 1, :, Array(binary_mask)]) ==
                  (val, val)

        end
    end
end

"""
    check_ocean_values_Y(Y, binary_mask; val = 0.0)

This function tests that every field stored in `Y` has all of
its values (where binary_mask == 1) equal to  `val`. Note that
this is meant to be used with the full land model (canopy,
snow, soil, soilco2).

Useful for checking if land model functions are updating the
values over the ocean.
"""
function check_ocean_values_Y(Y, binary_mask; val = 0.0)
    @test extrema(Array(parent(Y.soil.ϑ_l))[:, 1, 1, 1, Array(binary_mask)]) ==
          (val, val)
    @test extrema(Array(parent(Y.soil.θ_i))[:, 1, 1, 1, Array(binary_mask)]) ==
          (val, val)
    @test extrema(
        Array(parent(Y.soil.ρe_int))[:, 1, 1, 1, Array(binary_mask)],
    ) == (val, val)
    @test extrema(Array(parent(Y.snow.U))[1, 1, 1, Array(binary_mask)]) ==
          (val, val)
    @test extrema(Array(parent(Y.snow.S))[1, 1, 1, Array(binary_mask)]) ==
          (val, val)
    @test extrema(Array(parent(Y.snow.S_l))[1, 1, 1, Array(binary_mask)]) ==
          (val, val)
    @test extrema(
        Array(parent(Y.canopy.energy.T))[1, 1, 1, Array(binary_mask)],
    ) == (val, val)
    @test extrema(
        Array(parent(Y.canopy.hydraulics.ϑ_l.:1))[1, 1, 1, Array(binary_mask)],
    ) == (val, val)
end


context = ClimaComms.context()
nelements = (101, 15)
Δt = 450.0
FT = Float64
toml_dict = LP.create_toml_dict(FT)
earth_param_set = LP.LandParameters(toml_dict)

domain = ClimaLand.Domains.global_domain(
    FT;
    nelements = nelements,
    mask_threshold = FT(0.99),
);
surface_space = domain.space.surface;
start_date = DateTime(2008);
stop_date = start_date + Second(Δt)
forcing = ClimaLand.prescribed_forcing_era5(
    start_date,
    stop_date,
    domain.space.surface,
    toml_dict,
    FT;
    use_lowres_forcing = true,
    context,
)
LAI = ClimaLand.Canopy.prescribed_lai_modis(
    domain.space.surface,
    start_date,
    stop_date,
);

land = LandModel{FT}(forcing, LAI, toml_dict, domain, Δt);

@test domain == ClimaLand.get_domain(land)
@test ClimaComms.context(land) == ClimaComms.context()
@test ClimaComms.device(land) == ClimaComms.device()
@test ClimaLand.land_components(land) == (:soil, :snow, :soilco2, :canopy)

@testset "Total energy and water" begin
    Y, p, cds = initialize(land)
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
    area_index = p.canopy.biomass.area_index.leaf
    h_canopy = land.canopy.hydraulics.compartment_surfaces[end]
    ρ_ice = LP.ρ_cloud_ice(earth_param_set)
    ρ_liq = LP.ρ_cloud_liq(earth_param_set)
    int_cache = ClimaCore.Fields.zeros(domain.space.surface)
    ClimaCore.Operators.column_integral_definite!(
        int_cache,
        @. (ϑ_l0 + θ_i0 * ρ_ice / ρ_liq)
    )
    soil_exp = int_cache
    canopy_exp = ClimaCore.Fields.zeros(domain.space.surface)
    @. canopy_exp = area_index * h_canopy * ϑ0
    snow_exp = S0
    total_water = ClimaCore.Fields.zeros(domain.space.surface)
    cache = ClimaCore.Fields.zeros(domain.space.surface)
    ClimaLand.total_liq_water_vol_per_area!(total_water, land, Y, p, t0, cache)

    if pkgversion(ClimaCore) >= v"0.14.30"
        oceans = .~Array(parent(domain.space.surface.grid.mask.is_active))[:]
        continents = Array(parent(domain.space.surface.grid.mask.is_active))[:]
        expected = ClimaCore.Fields.zeros(domain.space.surface)
        @. expected = snow_exp + canopy_exp + soil_exp
        @test all(
            Array(parent(total_water))[1, 1, 1, continents] .≈
            Array(parent(expected))[1, 1, 1, continents],
        )
        @test all(
            Array(parent(total_water))[1, 1, 1, oceans] .≈
            Array(parent(expected))[1, 1, 1, oceans],
        )

        int_cache .*= 0
        ClimaCore.Operators.column_integral_definite!(int_cache, ρe_int0)
        soil_exp = int_cache
        @. canopy_exp =
            area_index * land.canopy.energy.parameters.ac_canopy * CTemp0
        snow_exp = U0
        total_energy = ClimaCore.Fields.zeros(domain.space.surface)
        ClimaLand.total_energy_per_area!(total_energy, land, Y, p, t0, cache)
        @. expected = snow_exp + canopy_exp + soil_exp
        @test all(
            Array(parent(total_energy))[1, 1, 1, continents] .≈
            Array(parent(expected))[1, 1, 1, continents],
        )
        @test all(
            Array(parent(total_energy))[1, 1, 1, oceans] .≈
            Array(parent(expected))[1, 1, 1, oceans],
        )
    else
        @test all(
            parent(total_water) .≈ parent(snow_exp .+ canopy_exp .+ soil_exp),
        )

        int_cache .*= 0
        ClimaCore.Operators.column_integral_definite!(int_cache, ρe_int0)
        soil_exp = int_cache
        canopy_exp =
            @. area_index * land.canopy.energy.parameters.ac_canopy * CTemp0
        snow_exp = U0
        total_energy = ClimaCore.Fields.zeros(domain.space.surface)
        ClimaLand.total_energy_per_area!(total_energy, land, Y, p, t0, cache)
        @test all(
            parent(total_energy) .≈ parent(snow_exp .+ canopy_exp .+ soil_exp),
        )
    end
end

if pkgversion(ClimaCore) >= v"0.14.30"
    @testset "Column integral mask awareness" begin
        Y, p, cds = initialize(land)
        Y.soil.ϑ_l .= land.soil.parameters.ν .+ FT(1e-3)
        fill!(parent(p.soil.is_saturated), FT(0.5)) # integrand
        @test extrema(p.soil.h∇) == (0.0, 0.0) # integral (0,0)
        @. p.soil.is_saturated = ClimaLand.Soil.is_saturated(
            Y.soil.ϑ_l + Y.soil.θ_i,
            land.soil.parameters.ν,
        )
        ClimaCore.Operators.column_integral_definite!(
            p.soil.h∇,
            p.soil.is_saturated,
        )
        @test maximum(p.soil.h∇) ≈ FT(50 * 1) # computed the integral over land ∫is_sat dz = 1 x ∫dz = 1 x 50m
        @test minimum(p.soil.h∇) ≈ FT(0.0) # did not compute an integral over the ocean, does not update
    end


    @testset "Mask of full land" begin
        Y, p, cds = initialize(land)
        Y .= 0
        surface_space = axes(Y.snow.U)
        subsurface_space = axes(Y.soil.ϑ_l)
        binary_mask = .~parent(surface_space.grid.mask.is_active)[:]
        # Test that the cache is zero over the ocean
        @info("testing initial cache")
        check_ocean_values_p(p, binary_mask)

        # Set initial conditions
        ic_path = ClimaLand.Artifacts.soil_ic_2008_50m_path(; context)
        set_initial_state! =
            ClimaLand.Simulations.make_set_initial_state_from_file(
                ic_path,
                land,
            )
        t0 = 0.0
        set_initial_state!(Y, p, t0, land)
        # Now, set the cache with physical values and make sure there are no NaNs, or values set over the ocean
        set_initial_cache! = make_set_initial_cache(land)
        set_initial_cache!(p, Y, t0)
        @info("testing set cache")
        check_ocean_values_p(p, binary_mask)

        # Check tendency functions do not update the state over the ocean
        dY = similar(Y)
        # Implicit tendency
        @. dY = 0
        imp_tendency! = make_imp_tendency(land)
        imp_tendency!(dY, Y, p, t0)
        @info("testing implicit tendency")
        check_ocean_values_Y(dY, binary_mask)

        # Explicit tendency
        @. dY = 0
        exp_tendency! = make_exp_tendency(land)
        exp_tendency!(dY, Y, p, t0)
        @info("testing explicit tendency")
        check_ocean_values_Y(dY, binary_mask)


        # Jacobian checks
        @info("testing Jacobian updates")

        jacobian! = ClimaLand.make_jacobian(land)
        jac_prototype = ClimaLand.FieldMatrixWithSolver(Y)

        # Check that the jacobian update respects the mask
        jacobian!(jac_prototype, Y, p, Δt, t0)
        (; matrix) = jac_prototype
        ∂ϑres∂ϑ = matrix[@name(soil.ϑ_l), @name(soil.ϑ_l)]
        @test extrema(
            Array(parent(∂ϑres∂ϑ.entries.:1))[:, 1, 1, 1, Array(binary_mask)],
        ) == (0.0, 0.0)
        @test extrema(
            Array(parent(∂ϑres∂ϑ.entries.:2))[:, 1, 1, 1, Array(binary_mask)],
        ) == (0.0, 0.0)
        @test extrema(
            Array(parent(∂ϑres∂ϑ.entries.:3))[:, 1, 1, 1, Array(binary_mask)],
        ) == (0.0, 0.0)

        ∂ρeres∂ρe = matrix[@name(soil.ρe_int), @name(soil.ρe_int)]
        @test extrema(
            Array(parent(∂ρeres∂ρe.entries.:1))[:, 1, 1, 1, Array(binary_mask)],
        ) == (0.0, 0.0)
        @test extrema(
            Array(parent(∂ρeres∂ρe.entries.:2))[:, 1, 1, 1, Array(binary_mask)],
        ) == (0.0, 0.0)
        @test extrema(
            Array(parent(∂ρeres∂ρe.entries.:3))[:, 1, 1, 1, Array(binary_mask)],
        ) == (0.0, 0.0)
        ∂ρeres∂ϑ = matrix[@name(soil.ρe_int), @name(soil.ϑ_l)]
        @test extrema(
            Array(parent(∂ρeres∂ϑ.entries.:1))[:, 1, 1, 1, Array(binary_mask)],
        ) == (0.0, 0.0)
        @test extrema(
            Array(parent(∂ρeres∂ϑ.entries.:2))[:, 1, 1, 1, Array(binary_mask)],
        ) == (0.0, 0.0)
        @test extrema(
            Array(parent(∂ρeres∂ϑ.entries.:3))[:, 1, 1, 1, Array(binary_mask)],
        ) == (0.0, 0.0)

        ∂Tres∂T = matrix[@name(canopy.energy.T), @name(canopy.energy.T)]
        @test extrema(
            Array(parent(∂Tres∂T.entries.:1))[1, 1, 1, Array(binary_mask)],
        ) == (0.0, 0.0)


        # Now carry out a solve of Jx = b, with x = Y, and b = 1
        b = similar(Y)
        fill!(parent(b), 1)
        x = deepcopy(Y)
        fill!(parent(x), 1)
        @test axes(x.soil.ϑ_l).grid.horizontal_grid.mask ==
              surface_space.grid.mask
        @test axes(b.soil.ϑ_l).grid.horizontal_grid.mask ==
              surface_space.grid.mask
        check_ocean_values_Y(x, binary_mask; val = 1.0)
        check_ocean_values_Y(b, binary_mask; val = 1.0)


        MatrixFields.field_matrix_solve!(
            jac_prototype.solver,
            x,
            jac_prototype.matrix,
            b,
        )
        check_ocean_values_Y(x, binary_mask; val = 1.0)

        # Take a step
        jac_kwargs = (; jac_prototype = jac_prototype, Wfact = jacobian!)
        prob = SciMLBase.ODEProblem(
            CTS.ClimaODEFunction(
                T_exp! = exp_tendency!,
                T_imp! = SciMLBase.ODEFunction(imp_tendency!; jac_kwargs...),
                dss! = ClimaLand.dss!,
            ),
            Y,
            (t0, t0 + Δt),
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
            saveat = [t0, t0 + Δt],
        )
        u = sol.u[end]
        check_ocean_values_Y(u, binary_mask)
    end
end

using ClimaComms
ClimaComms.@import_required_backends
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
using ClimaUtilities.TimeManager: ITime
using ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
import ClimaTimeSteppers
FT = Float64
context = ClimaComms.context()
toml_dict = LP.create_toml_dict(FT)
longlat = FT.((-118.68524, 37.06156));
domain =
    ClimaLand.Domains.Column(; zlim = FT.((-15, 0)), nelements = 15, longlat)
surface_space = domain.space.surface

start_date = DateTime(2008)
stop_date = start_date + Second(450)

atmos, radiation = ClimaLand.prescribed_forcing_era5(
    start_date,
    stop_date,
    surface_space,
    toml_dict,
    FT;
    use_lowres_forcing = true,
    context,
);
forcing = (; atmos, radiation);

@testset "No lake LandModel" begin
    LAI = TimeVaryingInput((t) -> 1.0)
    land = LandModel{FT}(
        forcing,
        LAI,
        toml_dict,
        domain,
        450.0;
        prognostic_land_components = (:canopy, :snow, :soil),
        conservation = true,
    )

    Y, p, cds = initialize(land)
    # Test that the variables are added
    @test :∫F_vol_e_dt ∈ propertynames(Y)
    @test :∫F_vol_liq_water_dt ∈ propertynames(Y)
    @test :total_water ∈ propertynames(p)
    @test :total_energy ∈ propertynames(p)

    t0 = ITime(0, Second(1), start_date)
    tf = ITime(450, Second(1), start_date)
    Δt = tf - t0
    set_ic! =
        ClimaLand.Simulations.make_set_initial_state_from_atmos_and_parameters(
            land,
        )

    set_ic!(Y, p, t0, land)
    set_initial_cache! = ClimaLand.make_set_initial_cache(land)
    set_initial_cache!(p, Y, t0)
    # Test that the cache is set correctly
    sfc_cache = p.scratch1 .* 0
    energy_cache = similar(sfc_cache) .* 0
    water_cache = similar(sfc_cache) .* 0
    ClimaLand.total_liq_water_vol_per_area!(
        water_cache,
        land,
        Y,
        p,
        t0,
        sfc_cache,
    )
    @test water_cache == p.total_water
    ClimaLand.total_energy_per_area!(energy_cache, land, Y, p, t0, sfc_cache)
    @test energy_cache == p.total_energy
    exp_tendency! = ClimaLand.make_exp_tendency(land)
    imp_tendency! = ClimaLand.make_compute_imp_tendency(land)
    update_jacobian! = ClimaLand.make_compute_jacobian(land)
    jacobian = ClimaLand.initialize_jacobian(Y)
    jac_kwargs = (; jac_prototype = jacobian, Wfact = update_jacobian!)
    T_imp! = ClimaTimeSteppers.ODEFunction(imp_tendency!; jac_kwargs...)
    CL_cache_imp! = ClimaLand.make_update_implicit_cache(land)
    cache_imp! = (Y, p, t) -> CL_cache_imp!(p, Y, t)
    problem = ClimaTimeSteppers.ODEProblem(
        ClimaTimeSteppers.ClimaODEFunction(
            T_exp! = exp_tendency!,
            T_imp! = T_imp!,
            dss! = ClimaLand.dss!,
            cache_imp! = cache_imp!,
        ),
        Y,
        (t0, tf),
        p,
    )
    timestepper = ClimaTimeSteppers.IMEXAlgorithm(
        ClimaTimeSteppers.ARS111(),
        ClimaTimeSteppers.NewtonsMethod(
            max_iters = 3,
            update_j = ClimaTimeSteppers.UpdateEvery(
                ClimaTimeSteppers.NewNewtonIteration,
            ),
        ),
    )
    _integrator =
        ClimaTimeSteppers.init(problem, timestepper; dt = Δt, adaptive = false)
    ClimaTimeSteppers.step!(_integrator)
    Y = deepcopy(_integrator.u)
    p = deepcopy(_integrator.p)
    ClimaLand.total_liq_water_vol_per_area!(
        p.total_water,
        land,
        Y,
        p,
        t0,
        sfc_cache,
    )
    ClimaLand.total_energy_per_area!(p.total_energy, land, Y, p, t0, sfc_cache)
    Δ_water = p.total_water .- water_cache
    @test parent(Y.∫F_vol_liq_water_dt)[1] .≈ parent(Δ_water)[1]
    Δ_energy = p.total_energy .- energy_cache
    @test parent(Y.∫F_vol_e_dt)[1] .≈ parent(Δ_energy)[1]
end


@testset "Lake pixel Conservation" begin
    LAI = TimeVaryingInput((t) -> 0.0)
    land = LandModel{FT}(
        forcing,
        LAI,
        toml_dict,
        domain,
        450.0;
        prognostic_land_components = (:canopy, :lake, :snow, :soil),
        conservation = true,
    )
    # make it a lake!
    land.lake.inland_water_mask .= 1

    Y, p, cds = initialize(land)
    t0 = ITime(0, Second(1), start_date)
    tf = ITime(450, Second(1), start_date)
    Δt = tf - t0
    set_ic! =
        ClimaLand.Simulations.make_set_initial_state_from_atmos_and_parameters(
            land,
        )

    set_ic!(Y, p, t0, land)
    set_initial_cache! = ClimaLand.make_set_initial_cache(land)
    set_initial_cache!(p, Y, t0)
    # Test that the cache is set correctly
    sfc_cache = p.scratch1 .* 0
    energy_cache = similar(sfc_cache) .* 0
    water_cache = similar(sfc_cache) .* 0
    ClimaLand.total_liq_water_vol_per_area!(
        water_cache,
        land,
        Y,
        p,
        t0,
        sfc_cache,
    )
    @test water_cache == p.total_water
    ClimaLand.total_energy_per_area!(energy_cache, land, Y, p, t0, sfc_cache)
    @test energy_cache == p.total_energy
    exp_tendency! = ClimaLand.make_exp_tendency(land)
    imp_tendency! = ClimaLand.make_compute_imp_tendency(land)
    update_jacobian! = ClimaLand.make_compute_jacobian(land)
    jacobian = ClimaLand.initialize_jacobian(Y)
    jac_kwargs = (; jac_prototype = jacobian, Wfact = update_jacobian!)
    T_imp! = ClimaTimeSteppers.ODEFunction(imp_tendency!; jac_kwargs...)
    CL_cache_imp! = ClimaLand.make_update_implicit_cache(land)
    cache_imp! = (Y, p, t) -> CL_cache_imp!(p, Y, t)
    problem = ClimaTimeSteppers.ODEProblem(
        ClimaTimeSteppers.ClimaODEFunction(
            T_exp! = exp_tendency!,
            T_imp! = T_imp!,
            dss! = ClimaLand.dss!,
            cache_imp! = cache_imp!,
        ),
        Y,
        (t0, tf),
        p,
    )
    timestepper = ClimaTimeSteppers.IMEXAlgorithm(
        ClimaTimeSteppers.ARS111(),
        ClimaTimeSteppers.NewtonsMethod(
            max_iters = 3,
            update_j = ClimaTimeSteppers.UpdateEvery(
                ClimaTimeSteppers.NewNewtonIteration,
            ),
        ),
    )
    _integrator =
        ClimaTimeSteppers.init(problem, timestepper; dt = Δt, adaptive = false)
    ClimaTimeSteppers.step!(_integrator)
    Y = deepcopy(_integrator.u)
    p = deepcopy(_integrator.p)
    ClimaLand.total_liq_water_vol_per_area!(
        p.total_water,
        land,
        Y,
        p,
        t0,
        sfc_cache,
    )
    ClimaLand.total_energy_per_area!(p.total_energy, land, Y, p, t0, sfc_cache)
    Δ_water = p.total_water .- water_cache
    @test parent(Y.∫F_vol_liq_water_dt)[1] .≈ parent(Δ_water)[1]
    Δ_energy = p.total_energy .- energy_cache
    @test parent(Y.∫F_vol_e_dt)[1] .≈ parent(Δ_energy)[1]
end

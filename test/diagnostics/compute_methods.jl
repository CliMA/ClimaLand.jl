using Test
using Dates
using ClimaLand
using ClimaLand: Domains, Soil, Canopy, Snow, Bucket
using ClimaLand.Diagnostics
import ClimaLand.Parameters as LP
import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput

# Every diagnostic a model declares possible must be defined for it and its
# compute function must run on the model's initialized state, both when
# allocating the output and when writing into a preallocated field.
function check_all_diagnostics(model, Y, p, t)
    empty!(Diagnostics.ALL_DIAGNOSTICS)
    names = Diagnostics.get_possible_diagnostics(model)
    @test allunique(names)
    Diagnostics.define_diagnostics!(model, names)
    @test Set(keys(Diagnostics.ALL_DIAGNOSTICS)) == Set(names)
    for name in names
        diag = Diagnostics.get_diagnostic_variable(name)
        out = diag.compute!(nothing, Y, p, t)
        @test out isa ClimaCore.Fields.Field
        out .= 0
        diag.compute!(out, Y, p, t)
        @test out isa ClimaCore.Fields.Field
    end
    short_names = Diagnostics.get_short_diagnostics(model)
    @test all(name -> name in names, short_names)
    return names
end

function initialize_with(set_ic!, model, t0)
    Y, p, cds = initialize(model)
    set_ic!(Y, p, t0, model)
    set_initial_cache! = make_set_initial_cache(model)
    set_initial_cache!(p, Y, t0)
    return Y, p
end

# Parameter-based initial conditions only exist in ClimaLand for LandModel and
# BucketModel (`Simulations.make_set_initial_state_from_atmos_and_parameters`).
# The standalone models use the minimal helpers below; the values are
# arbitrary and only need to keep the initial cache update finite.
function set_soil_ic!(Y, soil, ::Type{FT}) where {FT}
    Y.soil.ϑ_l .= FT(0.24)
    Y.soil.θ_i .= FT(0.0)
    ρc_s = Soil.volumetric_heat_capacity.(
        Y.soil.ϑ_l,
        Y.soil.θ_i,
        soil.parameters.ρc_ds,
        soil.parameters.earth_param_set,
    )
    Y.soil.ρe_int .= Soil.volumetric_internal_energy.(
        Y.soil.θ_i,
        ρc_s,
        FT(290.15),
        soil.parameters.earth_param_set,
    )
end

function set_canopy_ic!(Y, p, canopy, ::Type{FT}) where {FT}
    Y.canopy.hydraulics.ϑ_l .= canopy.hydraulics.parameters.ν
    Y.canopy.energy.T = FT(297.5)
    p.canopy.biomass.area_index.leaf .= FT(0.3)
    p.canopy.biomass.area_index.stem .= FT(0)
    p.canopy.biomass.area_index.root .= FT(0.3)
end

function set_soilco2_ic!(Y, ::Type{FT}) where {FT}
    Y.soilco2.CO2 .= FT(6e-5)
    Y.soilco2.O2 .= FT(0.08)
    Y.soilco2.SOC .= FT(5)
end

function set_snow_ic!(Y, ::Type{FT}) where {FT}
    Y.snow.S .= FT(0.01)
    Y.snow.S_l .= FT(0)
    Y.snow.U .= FT(-1e5)
end

FT = Float32
toml_dict = LP.create_toml_dict(FT)
domain = Domains.Column(;
    zlim = FT.((-1.0, 0.0)),
    nelements = 10,
    longlat = FT.((-118.1, 34.1)),
)
surface_domain = Domains.obtain_surface_domain(domain)
start_date = DateTime(2008)
stop_date = start_date + Day(3)
atmos, radiation = ClimaLand.prescribed_forcing_era5(
    start_date,
    stop_date,
    domain.space.surface,
    toml_dict,
    FT;
    use_lowres_forcing = true,
)
forcing = (; atmos, radiation)
dt = FT(900)
LAI = TimeVaryingInput((t) -> FT(1.0))
t0 = 0.0

@testset "EnergyHydrology diagnostics compute, FT = $FT" begin
    model = Soil.EnergyHydrology{FT}(domain, forcing, toml_dict)
    Y, p = initialize_with(model, t0) do Y, p, t0, model
        set_soil_ic!(Y, model, FT)
    end
    names = check_all_diagnostics(model, Y, p, t0)
    @test "swc" in names
    # A diagnostic that is not defined for a model raises an informative error
    @test_throws ErrorException Diagnostics.compute_snow_water_equivalent!(
        nothing,
        Y,
        p,
        t0,
        model,
    )
end

@testset "SoilCO2Model diagnostics compute, FT = $FT" begin
    T_soil = (z, t) -> eltype(z)(303)
    θ_l = (z, t) -> eltype(z)(0.3)
    hcm = Soil.vanGenuchten{FT}(; α = FT(0.1), n = FT(2))
    prescribed_met = Soil.Biogeochemistry.PrescribedMet{FT}(
        T_soil,
        θ_l,
        FT(0.6),
        FT(0.0),
        hcm,
    )
    drivers = Soil.Biogeochemistry.SoilDrivers(prescribed_met, atmos)
    model = Soil.Biogeochemistry.SoilCO2Model{FT}(domain, drivers, toml_dict)
    Y, p = initialize_with(model, t0) do Y, p, t0, model
        set_soilco2_ic!(Y, FT)
    end
    check_all_diagnostics(model, Y, p, t0)
end

@testset "CanopyModel diagnostics compute, FT = $FT" begin
    ground = ClimaLand.PrescribedGroundConditions{FT}()
    model = Canopy.CanopyModel{FT}(
        surface_domain,
        (; atmos, radiation, ground),
        LAI,
        toml_dict,
    )
    Y, p = initialize_with(model, t0) do Y, p, t0, model
        set_canopy_ic!(Y, p, model, FT)
    end
    check_all_diagnostics(model, Y, p, t0)
end

@testset "SnowModel diagnostics compute, FT = $FT" begin
    model = Snow.SnowModel(FT, surface_domain, forcing, toml_dict, dt)
    Y, p = initialize_with(model, t0) do Y, p, t0, model
        set_snow_ic!(Y, FT)
    end
    check_all_diagnostics(model, Y, p, t0)
end

@testset "BucketModel diagnostics compute, FT = $FT" begin
    bucket_atmos, bucket_rad =
        ClimaLand.prescribed_analytic_forcing(FT; toml_dict)
    albedo = Bucket.PrescribedBaregroundAlbedo{FT}(
        FT(0.8),
        (coordinate_point) -> 0.2,
        domain.space.surface,
    )
    parameters = Bucket.BucketModelParameters(
        toml_dict;
        albedo,
        z_0m = FT(1e-2),
        z_0b = FT(1e-3),
        τc = FT(1),
    )
    model = Bucket.BucketModel(;
        parameters,
        domain,
        atmosphere = bucket_atmos,
        radiation = bucket_rad,
    )
    set_ic! =
        ClimaLand.Simulations.make_set_initial_state_from_atmos_and_parameters(
            model,
        )
    Y, p = initialize_with(set_ic!, model, t0)
    check_all_diagnostics(model, Y, p, t0)
end

@testset "SoilCanopyModel diagnostics compute, FT = $FT" begin
    ground = ClimaLand.PrognosticGroundConditions{FT}()
    model = ClimaLand.SoilCanopyModel{FT}(
        (; atmos, radiation, ground),
        LAI,
        toml_dict,
        domain,
    )
    Y, p = initialize_with(model, t0) do Y, p, t0, model
        set_soil_ic!(Y, model.soil, FT)
        set_soilco2_ic!(Y, FT)
        set_canopy_ic!(Y, p, model.canopy, FT)
    end
    check_all_diagnostics(model, Y, p, t0)
end

@testset "LandModel diagnostics compute, FT = $FT" begin
    model = ClimaLand.LandModel{FT}(
        forcing,
        LAI,
        toml_dict,
        domain,
        dt;
        prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
    )
    set_ic! =
        ClimaLand.Simulations.make_set_initial_state_from_atmos_and_parameters(
            model,
        )
    Y, p = initialize_with(set_ic!, model, t0)
    names = check_all_diagnostics(model, Y, p, t0)
    @test "nee" in names
end

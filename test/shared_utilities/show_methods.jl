using Test
using Dates
using ClimaLand
using ClimaLand: Domains, Soil, Canopy, Bucket, Snow
using ClimaLand.Domains:
    Point, Column, Plane, HybridBox, SphericalShell, SphericalSurface
using ClimaLand.Soil:
    RichardsParameters,
    EnergyHydrologyParameters,
    WaterFluxBC,
    HeatFluxBC,
    WaterHeatBC,
    vanGenuchten
using ClimaLand.Soil.Runoff:
    NoRunoff, SurfaceRunoff, TOPMODELRunoff, TOPMODELSubsurfaceRunoff
import ClimaLand.Parameters as LP
import ClimaCore
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput

# Every custom show method should produce a multi-line text/plain form that
# names the type, a one-line compact form, and `summary` should match the
# compact form.
function check_show(x; typename = string(nameof(typeof(x))))
    out = sprint(show, MIME("text/plain"), x)
    @test occursin(typename, out)
    @test count(==('\n'), out) <= 12

    out2 = sprint(show, x)
    @test occursin(typename, out2)
    @test !occursin('\n', out2)
    out3 = sprint(show, MIME("text/plain"), x; context = :compact => true)
    @test out2 == out3

    out_summary = sprint(summary, x)
    @test occursin(typename, out_summary)
    @test !occursin('\n', out_summary)
    return out
end

struct DummyDomain{FT, S, F} <: Domains.AbstractDomain{FT}
    n::Int
    space::S
    fields::F
end

FT = Float32
toml_dict = LP.create_toml_dict(FT)
earth_param_set = LP.LandParameters(toml_dict)

@testset "Domain show methods, FT = $FT" begin
    column = Column(;
        zlim = FT.((-1.0, 0.0)),
        nelements = 5,
        dz_tuple = FT.((0.5, 0.05)),
    )
    out = check_show(column)
    @test occursin("mesh stretching", out)
    @test occursin("boundary names", out)

    plane = Plane(;
        xlim = FT.((0.0, 1.0)),
        ylim = FT.((0.0, 1.0)),
        nelements = (2, 2),
        longlat = FT.((0.0, 0.0)),
    )
    out = check_show(plane)
    @test occursin("longlat center", out)

    box = HybridBox(;
        xlim = FT.((0.0, 1.0)),
        ylim = FT.((0.0, 1.0)),
        zlim = FT.((-1.0, 0.0)),
        nelements = (2, 2, 4),
        longlat = FT.((0.0, 0.0)),
        dz_tuple = FT.((0.5, 0.05)),
    )
    out = check_show(box)
    @test occursin("longlat center", out)
    @test occursin("mesh stretching", out)

    shell = SphericalShell(;
        radius = FT(100),
        depth = FT(30),
        nelements = (2, 4),
        dz_tuple = FT.((10.0, 1.0)),
    )
    out = check_show(shell)
    @test occursin("mesh stretching", out)

    surface = SphericalSurface(; radius = FT(100), nelements = 2)
    out = check_show(surface)
    @test occursin("radius", out)

    # A domain without a dedicated show method falls back to the generic
    # AbstractDomain method, which lists every field except the spaces.
    dummy = DummyDomain{FT, typeof(column.space), typeof(column.fields)}(
        3,
        column.space,
        column.fields,
    )
    out = check_show(dummy)
    @test occursin("n: 3", out)
    @test !occursin("space", out)
    out2 = sprint(show, dummy)
    @test occursin("n=3", out2)
end

@testset "Parameter show methods, FT = $FT" begin
    out = check_show(earth_param_set)
    @test occursin("physical constants", out)

    hcm = vanGenuchten{FT}(; α = FT(2.6), n = FT(2.0))
    richards = RichardsParameters(;
        ν = FT(0.495),
        hydrology_cm = hcm,
        K_sat = FT(1e-5),
        S_s = FT(1e-3),
        θ_r = FT(0),
    )
    out = check_show(richards)
    @test occursin("vanGenuchten", out)
    @test occursin("0.495", out)

    energy_hydrology = EnergyHydrologyParameters(
        toml_dict;
        ν = FT(0.495),
        ν_ss_om = FT(0.0),
        ν_ss_quartz = FT(1.0),
        ν_ss_gravel = FT(0.0),
        hydrology_cm = hcm,
        K_sat = FT(1e-5),
        S_s = FT(1e-3),
        θ_r = FT(0),
    )
    out = check_show(energy_hydrology)
    @test occursin("vanGenuchten", out)
    @test occursin("albedo", out)

    # Spatially varying parameters print the type name, not the values
    column = Column(; zlim = FT.((-1.0, 0.0)), nelements = 5)
    ν_field = ClimaCore.Fields.zeros(column.space.subsurface) .+ FT(0.4)
    hcm_field = map(_ -> hcm, ν_field)
    richards_field = RichardsParameters(;
        ν = ν_field,
        hydrology_cm = hcm_field,
        K_sat = ν_field .* FT(1e-5),
        S_s = ν_field .* FT(1e-3),
        θ_r = ν_field .* FT(0),
    )
    out = check_show(richards_field)
    @test occursin("Field", out)
    @test occursin("vanGenuchten", out)

    snow = Snow.SnowParameters(toml_dict, FT(450))
    out = check_show(snow)
    @test occursin("timestep", out)

    bucket_domain =
        SphericalShell(; radius = FT(100), depth = FT(3.5), nelements = (1, 10))
    albedo = Bucket.PrescribedBaregroundAlbedo{FT}(
        FT(0.8),
        (coordinate_point) -> 0.2,
        bucket_domain.space.surface,
    )
    out = check_show(albedo)
    @test occursin("α_snow", out)
    bucket = Bucket.BucketModelParameters(
        toml_dict;
        albedo,
        z_0m = FT(1e-2),
        z_0b = FT(1e-3),
        τc = FT(1.0),
    )
    out = check_show(bucket)
    @test occursin("bucket capacity", out)
    @test occursin("PrescribedBaregroundAlbedo", out)

    point = Point(; z_sfc = FT(0), longlat = FT.((-180, 1)))
    radiation_parameters = (;
        α_PAR_leaf = FT(0.1),
        α_NIR_leaf = FT(0.4),
        Ω = FT(1),
        G_Function = Canopy.ConstantGFunction(FT(0.5)),
    )
    two_stream = Canopy.TwoStreamModel{FT}(point, toml_dict)
    out = check_show(two_stream.parameters)
    @test occursin("clumping index", out)
    @test occursin("CLMGFunction", out)
end

@testset "Boundary condition and runoff show methods, FT = $FT" begin
    water_bc = WaterFluxBC((p, t) -> 0.0)
    out = check_show(water_bc)
    @test occursin("<prescribed function>", out)

    heat_bc = HeatFluxBC((p, t) -> 0.0)
    water_heat_bc = WaterHeatBC(; water = water_bc, heat = heat_bc)
    out = check_show(water_heat_bc)
    @test occursin("WaterFluxBC", out)
    @test occursin("HeatFluxBC", out)

    check_show(NoRunoff())
    check_show(SurfaceRunoff())

    column = Column(; zlim = FT.((-1.0, 0.0)), nelements = 5)
    f_max = ClimaCore.Fields.zeros(column.space.surface) .+ FT(0.5)
    topmodel = TOPMODELRunoff{FT}(; f_over = FT(3.28), f_max, R_sb = FT(1e-4))
    out = check_show(topmodel)
    @test occursin("f_max: Field", out)
    @test occursin("TOPMODELSubsurfaceRunoff", out)

    subsurface = TOPMODELSubsurfaceRunoff{FT}(FT(1e-4), FT(3.28))
    out = check_show(subsurface)
    @test occursin("R_sb", out)
    @test occursin("explicit=false", out)
end

@testset "Driver show methods, FT = $FT" begin
    start_date = DateTime(2005)
    f = TimeVaryingInput((t) -> 10.0)
    g = TimeVaryingInput((t) -> 20.0)
    atmos = ClimaLand.PrescribedAtmosphere(
        f,
        f,
        f,
        f,
        f,
        f,
        start_date,
        FT(1),
        toml_dict;
        c_co2 = g,
    )
    out = check_show(atmos)
    @test occursin("reference height", out)
    @test occursin("CO2", out)

    radiation = ClimaLand.PrescribedRadiativeFluxes(
        FT,
        f,
        f,
        start_date;
        cosθs = (t, s) -> 1.0,
        frac_diff = g,
    )
    out = check_show(radiation)
    @test occursin("cosθs: prescribed", out)
    @test occursin("frac_diff: prescribed", out)

    radiation_default =
        ClimaLand.PrescribedRadiativeFluxes(FT, f, g, start_date; toml_dict)
    out = check_show(radiation_default)
    @test occursin("computed from date and location", out)
    @test occursin("computed empirically", out)
end

@testset "Canopy component show methods, FT = $FT" begin
    point = Point(; z_sfc = FT(0), longlat = FT.((-180, 1)))
    radiation_parameters = (;
        α_PAR_leaf = FT(0.1),
        α_NIR_leaf = FT(0.4),
        Ω = FT(1),
        G_Function = Canopy.ConstantGFunction(FT(0.5)),
    )
    for component in (
        Canopy.BeerLambertModel{FT}(point, toml_dict; radiation_parameters),
        Canopy.TwoStreamModel{FT}(point, toml_dict),
        Canopy.BigLeafEnergyModel{FT}(toml_dict),
        Canopy.AutotrophicRespirationModel{FT}(toml_dict),
    )
        out = check_show(component)
        @test occursin("parameters", out)
    end
end

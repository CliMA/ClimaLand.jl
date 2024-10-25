using Test
using ClimaLand
using ClimaLand.Diagnostics: @with_error

@test isdefined(ClimaLand.Diagnostics, :compute_sw_albedo!)

@test !hasmethod(
    ClimaLand.Diagnostics.compute_sw_albedo!,
    (Any, Any, Any, Any, Any),
)

# Define some diagnostics for a DummyModel

@test ClimaLand.Diagnostics.ALL_DIAGNOSTICS isa Dict
@test length(ClimaLand.Diagnostics.ALL_DIAGNOSTICS) == 0
struct DummyModel end
struct DummyModel2 end
ClimaLand.Diagnostics.@diagnostic_compute "sw_albedo" Union{
    DummyModel,
    DummyModel2,
} p.foo.bar

ClimaLand.Diagnostics.add_diagnostic_variable!(
    short_name = "alpha",
    long_name = "Albedo",
    standard_name = "albedo",
    units = "",
    compute! = (out, Y, p, t) -> compute_sw_albedo!(out, Y, p, t, land_model),
)

@test length(ClimaLand.Diagnostics.ALL_DIAGNOSTICS) == 1

ClimaLand.Diagnostics.define_diagnostics!(DummyModel())

# Just to trigger the error
out = Y = p = t = land_model = nothing

@test_throws ErrorException("Cannot compute albedo with model = Nothing") ClimaLand.Diagnostics.compute_albedo!(
    out,
    Y,
    p,
    t,
    land_model,
)

@test_throws ErrorException ClimaLand.Diagnostics.get_diagnostic_variable("Foo")

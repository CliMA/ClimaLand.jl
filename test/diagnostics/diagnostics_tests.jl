using Test
using ClimaLand

@test isdefined(ClimaLand.Diagnostics, :compute_albedo!)

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

import ClimaParams as CP
using Test

override_file = joinpath("toml", "tags.toml")
toml_dict = CP.create_toml_dict(Float32; override_file)

beljaars_param_names = sort([
    "coefficient_a_m_beljaars",
    "coefficient_b_h_beljaars",
    "coefficient_c_h_beljaars",
    "coefficient_c_m_beljaars",
    "coefficient_d_h_beljaars",
    "most_stability_parameter_beljaars",
    "coefficient_d_m_beljaars",
    "coefficient_b_m_beljaars",
    "prandtl_number_0_beljaars",
    "coefficient_a_h_beljaars",
    "most_stability_exponent_beljaars",
])

beljaars_params = CP.get_parameter_values(toml_dict, beljaars_param_names)

@testset "Tags" begin

    @test pairs(beljaars_params) ==
          pairs(CP.get_tagged_parameter_values(toml_dict, ["bel jaa*rs"]))

    @test beljaars_param_names ==
          sort(CP.get_tagged_parameter_names(toml_dict, ["bel jaa*rs"]))

end

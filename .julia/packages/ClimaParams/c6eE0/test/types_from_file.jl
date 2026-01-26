using Test

import ClimaParams as CP

path_to_params = joinpath(@__DIR__, "toml", "typed_parameters.toml")

@testset "parameter types from file" begin
    param_names = [
        "int_array_param",
        "int_param",
        "string_param",
        "string_array_param",
        "ft_array_param",
        "bool_param",
        "light_speed",
    ]
    for FT in (Float32, Float64)
        td = CP.create_toml_dict(FT; override_file = path_to_params)

        nt = CP.get_parameter_values(td, param_names)
        @test typeof(nt.string_param) == String
        @test eltype(nt.string_array_param) == String

        @test typeof(nt.int_param) == Int
        @test eltype(nt.int_array_param) == Int

        @test typeof(nt.light_speed) == FT
        @test eltype(nt.ft_array_param) == FT

        @test_broken CP.get_parameter_values(td, "untyped_param")

        @test_throws ErrorException CP.get_parameter_values(td, "badtype_param")
    end
end

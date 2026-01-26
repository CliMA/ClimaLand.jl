using Test

import ClimaParams as CP

@testset "Merge correctness" begin

    FT = Float32
    toml_1 = ("toml/merge1.toml")
    toml_2 = ("toml/merge2.toml")

    toml_dict_1 = CP.create_toml_dict(FT, default_file = toml_1)
    toml_dict_2 = CP.create_toml_dict(FT, default_file = toml_2)

    # Test individual TOMLs
    (; a) = CP.get_parameter_values(toml_dict_1, "a")
    @test a == 0.0
    @test a isa FT
    (; a) = CP.get_parameter_values(toml_dict_2, "a")
    @test a == 2
    @test a isa Int

    # Test direct access
    @test toml_dict_1["a"] == 0.0
    @test toml_dict_1["a"] isa FT
    @test toml_dict_2["a"] == 2
    @test toml_dict_2["a"] isa Int

    # Test if parameter is logged by getindex
    @test isnothing(CP.check_override_parameter_usage(toml_dict_1, true))
    @test isnothing(CP.check_override_parameter_usage(toml_dict_2, true))

    # Test merging
    merged_dict = CP.merge_toml_files([toml_1, toml_2]; override = true)
    merged_toml_dict = CP.create_toml_dict(FT, default_file = merged_dict)
    (; a, b) = CP.get_parameter_values(merged_toml_dict, ["a", "b"])
    @test a == 2
    @test a isa Int
    @test b == 3.0
    @test b isa FT

    # Swap merge order
    merged_dict = CP.merge_toml_files([toml_2, toml_1]; override = true)
    merged_toml_dict = CP.create_toml_dict(FT, default_file = merged_dict)
    (; a, b) = CP.get_parameter_values(merged_toml_dict, ["a", "b"])
    @test a == 0.0
    @test a isa FT
    @test b == 3.0
    @test b isa FT
end

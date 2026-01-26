module ParameterBoxes

using Test

import ClimaParams as CP

Base.@kwdef struct ParameterBox{FT}
    molmass_dryair::FT
    universal_gas_constant::FT
    new_parameter::FT
end

# Derived parameters
R_d(pb::ParameterBox) = pb.universal_gas_constant / pb.molmass_dryair

name_map = Dict(
    :universal_gas_constant => :universal_gas_constant,
    :molar_mass_dry_air => :molmass_dryair,
    :param_2 => :new_parameter,
)

@testset "Example use case: parameter box" begin

    # [1.] read from file
    toml_file = joinpath(@__DIR__, "toml", "parambox.toml")
    toml_dict = CP.create_toml_dict(Float64; override_file = toml_file)

    # [2.] build
    param_set = CP.create_parameter_struct(ParameterBox, toml_dict, name_map)

    # [3.] log & checks(with warning)
    mktempdir(@__DIR__) do path
        logfilepath = joinpath(path, "logfilepath.toml")
        @test_logs (:warn,) CP.log_parameter_information(toml_dict, logfilepath)
    end

    # [4.] use
    # from default
    @test param_set.molmass_dryair ≈ 0.02897
    # overridden default
    @test param_set.universal_gas_constant ≈ 4.0
    # derived in constructor
    @test R_d(param_set) ≈
          param_set.universal_gas_constant / param_set.molmass_dryair
    # from toml
    @test param_set.new_parameter ≈ 19.99

end

end # end module

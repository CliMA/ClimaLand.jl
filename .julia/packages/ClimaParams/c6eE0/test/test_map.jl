using Test
import ClimaParams as CP
import Thermodynamics.Parameters.ThermodynamicsParameters

@testset "name map tests" begin
    toml_dict = CP.create_toml_dict(Float64)

    @testset "Smoke test" begin
        symbol_dict = Dict(
            :gravitational_acceleration => :g,
            :angular_velocity_planet_rotation => :omega,
        )
        string_dict = Dict(
            "gravitational_acceleration" => "g",
            "angular_velocity_planet_rotation" => "omega",
        )
        symbol_pairs = (
            :gravitational_acceleration => :g,
            :angular_velocity_planet_rotation => :omega,
        )
        string_pairs = (
            "gravitational_acceleration" => "g",
            "angular_velocity_planet_rotation" => "omega",
        )
        for name_map in (
            symbol_dict,
            string_dict,
            symbol_pairs,
            string_pairs,
            [symbol_pairs...],
            [string_pairs...],
        )
            pairs = CP.get_parameter_values(toml_dict, name_map)
            @test pairs.g == 9.81
            @test pairs.omega == 7.2921159e-5
        end

        # Test varargs
        pairs = CP.get_parameter_values(
            toml_dict,
            :gravitational_acceleration => :g,
            :angular_velocity_planet_rotation => :omega,
        )

        @test pairs.g == 9.81
        @test pairs.omega == 7.2921159e-5

        # Test tuple
        pairs = CP.get_parameter_values(
            toml_dict,
            (
                :gravitational_acceleration => :g,
                :angular_velocity_planet_rotation => :omega,
            ),
        )
        @test pairs.g == 9.81
        @test pairs.omega == 7.2921159e-5
    end


    ThermodynamicsParameterMap = (;
        :temperature_min_reference => :T_min_ref,
        :entropy_water_vapor => :entropy_water_vapor,
        :entropy_dry_air => :entropy_dry_air,
        :entropy_reference_temperature => :entropy_reference_temperature,
        :temperature_saturation_adjustment_max => :T_max,
        :molar_mass_dry_air => :molmass_dryair,
        :pow_icenuc => :pow_icenuc,
        :temperature_triple_point => :T_triple,
        :adiabatic_exponent_dry_air => :kappa_d,
        :pressure_triple_point => :press_triple,
        :thermodynamics_temperature_reference => :T_0,
        :temperature_water_freeze => :T_freeze,
        :isobaric_specific_heat_ice => :cp_i,
        :latent_heat_sublimation_at_reference => :LH_s0,
        :isobaric_specific_heat_vapor => :cp_v,
        :molar_mass_water => :molmass_water,
        :mean_sea_level_pressure => :MSLP,
        :isobaric_specific_heat_liquid => :cp_l,
        :latent_heat_vaporization_at_reference => :LH_v0,
        :temperature_saturation_adjustment_min => :T_min,
        # :temperature_saturation_adjustment_init_min => :T_init_min, # will need updated
        :universal_gas_constant => :gas_constant,
        :temperature_surface_reference => :T_surf_ref,
        :gravitational_acceleration => :grav,
        :temperature_homogenous_nucleation => :T_icenuc,
        :potential_temperature_reference_pressure => :p_ref_theta,
    )

    # Example function for thermo params - essentially a specific constuctor of `create_parameter_struct`
    function thermo_params(toml_dict)
        params = CP.get_parameter_values(toml_dict, ThermodynamicsParameterMap)
        FT = CP.float_type(toml_dict)
        ThermodynamicsParameters{FT}(; params...)
    end

    @testset "Test `create_parameter_struct`" begin
        _thermo_params = thermo_params(toml_dict)
        @test _thermo_params == CP.create_parameter_struct(
            ThermodynamicsParameters,
            toml_dict,
            ThermodynamicsParameterMap,
        )
    end

end

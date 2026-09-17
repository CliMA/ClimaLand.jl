using Test
using Dates
using ClimaLand
using ClimaLand: Domains
import ClimaLand.Parameters as LP
import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore
import ClimaUtilities.TimeVaryingInputs:
    LinearInterpolation, PeriodicCalendar, evaluate!
import Thermodynamics

FT = Float32
toml_dict = LP.create_toml_dict(FT)
earth_param_set = LP.LandParameters(toml_dict)
thermo_params = LP.thermodynamic_parameters(earth_param_set)

@testset "Perturbed specific humidity from dewpoint, FT = $FT" begin
    T_air = FT(290)
    T_dew = FT(285)
    P_air = FT(101325)
    e_sat(T) = Thermodynamics.saturation_vapor_pressure(
        thermo_params,
        T,
        Thermodynamics.Liquid(),
    )
    rh = e_sat(T_dew) / e_sat(T_air)
    q_expected = Thermodynamics.q_vap_from_RH(
        thermo_params,
        P_air,
        T_air,
        rh,
        Thermodynamics.Liquid(),
    )

    q0 = ClimaLand.perturbed_temp_specific_humidity_from_dewpoint(
        T_dew,
        T_air,
        P_air,
        earth_param_set,
        FT(0),
    )
    @test q0 ≈ q_expected
    # Warmer air at fixed relative humidity holds more water
    q_warm = ClimaLand.perturbed_temp_specific_humidity_from_dewpoint(
        T_dew,
        T_air,
        P_air,
        earth_param_set,
        FT(2),
    )
    @test q_warm > q0
    # Data of a different float type are converted to the simulation type
    q64 = ClimaLand.perturbed_temp_specific_humidity_from_dewpoint(
        285.0,
        290.0,
        101325.0,
        earth_param_set,
        0.0,
    )
    @test q64 isa FT
    @test q64 ≈ q0

    q_rh0 = ClimaLand.perturbed_rh_specific_humidity_from_dewpoint(
        T_dew,
        T_air,
        P_air,
        earth_param_set,
        FT(0),
    )
    @test q_rh0 ≈ q0
    q_rh_up = ClimaLand.perturbed_rh_specific_humidity_from_dewpoint(
        T_dew,
        T_air,
        P_air,
        earth_param_set,
        FT(0.1),
    )
    @test q_rh_up > q_rh0
    # The perturbed relative humidity is clipped to (0, 1]
    q_rh_sat = ClimaLand.perturbed_rh_specific_humidity_from_dewpoint(
        T_dew,
        T_air,
        P_air,
        earth_param_set,
        FT(1),
    )
    @test q_rh_sat ≈ Thermodynamics.q_vap_from_RH(
        thermo_params,
        P_air,
        T_air,
        FT(1),
        Thermodynamics.Liquid(),
    )
    q_rh_dry = ClimaLand.perturbed_rh_specific_humidity_from_dewpoint(
        T_dew,
        T_air,
        P_air,
        earth_param_set,
        FT(-1),
    )
    @test 0 < q_rh_dry < q_rh0
end

@testset "Perturbed ERA5 forcing, FT = $FT" begin
    era5_path = ClimaLand.Artifacts.era5_land_forcing_data2008_lowres_path()
    domain = Domains.Column(;
        zlim = FT.((-1.0, 0.0)),
        nelements = 4,
        longlat = FT.((-118.1, 34.1)),
    )
    surface_space = domain.space.surface
    start_date = DateTime(2008)
    t = 0.0
    function evaluate(input)
        field = ClimaCore.Fields.zeros(surface_space)
        evaluate!(field, input, t)
        return field
    end

    ΔT = FT(2)
    reference = ClimaLand.prescribed_perturbed_temperature_era5(
        era5_path,
        surface_space,
        start_date,
        toml_dict,
        FT(0),
        FT,
    )
    warmer = ClimaLand.prescribed_perturbed_temperature_era5(
        era5_path,
        surface_space,
        start_date,
        toml_dict,
        ΔT,
        FT;
        max_wind_speed = FT(1),
        time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
    )
    @test reference.atmos isa ClimaLand.PrescribedAtmosphere
    @test reference.radiation isa ClimaLand.PrescribedRadiativeFluxes

    T_ref = evaluate(reference.atmos.T)
    T_warm = evaluate(warmer.atmos.T)
    @test all(parent(T_warm) .≈ parent(T_ref) .+ ΔT)
    # Fixed relative humidity: warmer air holds more water
    @test all(
        parent(evaluate(warmer.atmos.q)) .> parent(evaluate(reference.atmos.q)),
    )
    # Downwelling longwave shifts by the climate sensitivity a = 2 W/m^2/K
    LW_ref = evaluate(reference.radiation.LW_d)
    LW_warm = evaluate(warmer.radiation.LW_d)
    @test all(parent(LW_warm) .≈ parent(LW_ref) .+ 2ΔT)
    @test all(parent(evaluate(warmer.atmos.u)) .<= FT(1))
    frac_diff = evaluate(reference.radiation.frac_diff)
    @test all(FT(0) .<= parent(frac_diff) .<= FT(1))
    @test all(parent(evaluate(reference.atmos.liquid_precip)) .<= FT(0))
    @test all(parent(evaluate(reference.atmos.P)) .> FT(0))

    Δrh = FT(0.1)
    reference_rh = ClimaLand.prescribed_perturbed_rh_era5(
        era5_path,
        surface_space,
        start_date,
        toml_dict,
        FT(0),
        FT,
    )
    moister = ClimaLand.prescribed_perturbed_rh_era5(
        era5_path,
        surface_space,
        start_date,
        toml_dict,
        Δrh,
        FT;
        max_wind_speed = FT(1),
        time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
    )
    # The relative humidity perturbation leaves temperature and radiation unchanged
    @test parent(evaluate(moister.atmos.T)) ==
          parent(evaluate(reference_rh.atmos.T))
    @test parent(evaluate(reference_rh.atmos.T)) == parent(T_ref)
    @test parent(evaluate(moister.radiation.LW_d)) == parent(LW_ref)
    @test all(
        parent(evaluate(moister.atmos.q)) .>
        parent(evaluate(reference_rh.atmos.q)),
    )
    @test parent(evaluate(reference_rh.atmos.q)) ≈
          parent(evaluate(reference.atmos.q))
end

using Test
using Dates
using StaticArrays
using Insolation
import ClimaComms
ClimaComms.@import_required_backends
using ClimaCore
import ClimaParams as CP
using ClimaLand
using ClimaLand.Canopy
using ClimaLand.Canopy.PlantHydraulics
using ClimaLand.Domains: Point, Plane

import ClimaLand
import ClimaLand.Parameters as LP

for FT in (Float32, Float64)
    # Domain information
    cmax = FT(10)
    cmin = FT(0)
    nelems = 10
    longlat = (FT(-118.14), FT(34.15))
    pt = Point(; z_sfc = FT(0), longlat)
    plane = Plane(;
        xlim = (cmin, cmax),
        ylim = (cmin, cmax),
        nelements = (nelems, nelems),
        longlat,
    )
    domains = [pt, plane]

    # Parameters
    toml_dict = LP.create_toml_dict(FT)
    earth_param_set = LP.LandParameters(toml_dict)

    # Time information
    t0 = 0.0
    start_date = DateTime(2005)
    Δt = FT(180.0)

    # Radiation forcing
    SW_d = TimeVaryingInput((t) -> eltype(t)(20.0))
    LW_d = TimeVaryingInput((t) -> eltype(t)(20.0))
    zenith_angle =
        (t, s) -> default_zenith_angle(
            t,
            s;
            insol_params = earth_param_set.insol_params,
            latitude = FT(40.0),
            longitude = FT(-120.0),
        )

    radiation = PrescribedRadiativeFluxes(
        FT,
        SW_d,
        LW_d,
        start_date;
        θs = zenith_angle,
        toml_dict = toml_dict,
    )
    # Atmos forcing
    precip = TimeVaryingInput((t) -> eltype(t)(0)) # no precipitation
    T_atmos = TimeVaryingInput((t) -> eltype(t)(290.0))
    u_atmos = TimeVaryingInput((t) -> eltype(t)(2.0))
    q_atmos = TimeVaryingInput((t) -> eltype(t)(0.011))
    h_atmos = FT(3)
    P_atmos = TimeVaryingInput((t) -> eltype(t)(101325))
    atmos = ClimaLand.PrescribedAtmosphere(
        precip,
        precip,
        T_atmos,
        u_atmos,
        q_atmos,
        P_atmos,
        start_date,
        h_atmos,
        toml_dict,
    )
    # Ground forcing
    ground = PrescribedGroundConditions{FT}()

    # Hydraulics model
    LAI_value = FT(8)
    LAI = TimeVaryingInput(t -> LAI_value) # m2 [leaf] m-2 [ground]
    RAI = FT(1)
    SAI = FT(1)
    ai_parameterization =
        PlantHydraulics.PrescribedSiteAreaIndex{FT}(LAI, SAI, RAI)
    n_stem = Int64(2) # number of stem elements
    n_leaf = Int64(1) # number of leaf elements
    h_stem = h_leaf = FT(1)

    @testset "Canopy model total energy and water, FT = $FT" begin
        for domain in domains
            hydraulics = PlantHydraulics.PlantHydraulicsModel{FT}(
                domain,
                LAI,
                toml_dict;
                n_stem,
                n_leaf,
                h_stem,
                h_leaf,
                SAI,
                RAI,
                ai_parameterization,
            )
            canopy = ClimaLand.Canopy.CanopyModel{FT}(
                domain,
                (; radiation, atmos, ground),
                LAI,
                toml_dict;
                hydraulics,
            )

            Y, p, cds = initialize(canopy)
            ϑ0 = canopy.hydraulics.parameters.ν / 2
            Temp0 = FT(290.5)

            Y.canopy.hydraulics.ϑ_l.:1 .= ϑ0
            Y.canopy.hydraulics.ϑ_l.:2 .= ϑ0
            Y.canopy.hydraulics.ϑ_l.:3 .= ϑ0
            Y.canopy.energy.T .= Temp0

            set_initial_cache! = make_set_initial_cache(canopy)
            set_initial_cache!(p, Y, t0)
            total_energy = ClimaCore.Fields.zeros(domain.space.surface)
            ClimaLand.total_energy_per_area!(total_energy, canopy, Y, p, t0)
            @test all(
                parent(total_energy) .≈
                (LAI_value + SAI) * canopy.energy.parameters.ac_canopy * Temp0,
            )

            total_water = ClimaCore.Fields.zeros(domain.space.surface)
            ClimaLand.total_liq_water_vol_per_area!(
                total_water,
                canopy,
                Y,
                p,
                t0,
            )
            @test all(
                parent(total_water) .≈
                (LAI_value * n_leaf * h_leaf * ϑ0 + SAI * n_stem * h_stem * ϑ0),
            )
        end
    end
end

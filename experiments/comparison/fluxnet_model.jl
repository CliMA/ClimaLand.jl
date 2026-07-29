# The `run_fluxnet.jl` model, shared by the scripts in this directory so that a comparison and
# a benchmark provably step the same thing: an integrated `LandModel` of soil (with root
# extraction and surface runoff), soil CO2 biogeochemistry, a canopy (two-stream radiative
# transfer, Medlyn conductance, Farquhar photosynthesis, plant hydraulics, big-leaf energy,
# prescribed biomass with MODIS LAI), and snow, started from the site's observed initial
# conditions.
#
# `fluxnet_site_setup` collects everything that depends only on the site; `build_fluxnet_sim`
# turns a domain and a forcing into a `LandSimulation`, so two builds can differ in nothing
# else.

using Dates
using ClimaLand
using ClimaLand.Canopy
using ClimaLand.Snow
using ClimaLand.Soil
using ClimaLand.Soil.Biogeochemistry
import ClimaLand
import ClimaLand.Simulations: LandSimulation
import ClimaLand.FluxnetSimulations as FluxnetSimulations

"""
    fluxnet_site_setup(::Type{FT}, site_ID) where {FT}

Return a NamedTuple of everything about `site_ID` that a build shares with any other build of
the same site: the vertical mesh (`dz_tuple`, `nelements`, `zmin`, `zmax`), the location
(`time_offset`, `lat`, `long`), the flux tower height `atmos_h`, every soil and plant
parameter from `get_parameters`, the `start_date`/`stop_date` spanned by the site's data, the
location-based `maxLAI` with the `RAI` derived from it, and the `run_fluxnet.jl` timestep
`dt`.

Pass the result to [`build_fluxnet_sim`](@ref).
"""
function fluxnet_site_setup(::Type{FT}, site_ID) where {FT}
    site_val = Val(FluxnetSimulations.replace_hyphen(site_ID))
    (; dz_tuple, nelements, zmin, zmax) =
        FluxnetSimulations.get_domain_info(FT, site_val)
    (; time_offset, lat, long) = FluxnetSimulations.get_location(FT, site_val)
    (; atmos_h) = FluxnetSimulations.get_fluxtower_height(FT, site_val)
    parameters = FluxnetSimulations.get_parameters(FT, site_val)

    dt = FT(450) # 7.5 minutes
    (start_date, stop_date) =
        FluxnetSimulations.get_data_dates(site_ID, time_offset)
    maxLAI = FluxnetSimulations.get_maxLAI_at_site(start_date, lat, long)

    return (;
        site_ID,
        dz_tuple,
        nelements,
        zmin,
        zmax,
        time_offset,
        lat,
        long,
        atmos_h,
        dt,
        start_date,
        stop_date,
        maxLAI,
        RAI = maxLAI * parameters.f_root_to_shoot,
        parameters...,
    )
end

"""
    build_fluxnet_sim(land_domain, forcing, site, toml_dict;
                      user_callbacks = nothing,
                      diagnostics = nothing)

Build the `run_fluxnet.jl` `LandSimulation` on `land_domain`, driven by `forcing` (a
`(; atmos, radiation)` NamedTuple), for the site described by `site` from
[`fluxnet_site_setup`](@ref).

`land_domain` and `forcing` are the only things meant to vary between builds: the parameters,
vertical mesh, LAI, timestep, and initial conditions all come from `site`.

`user_callbacks` is forwarded to `LandSimulation` only when it is given, so the default
`NaNCheckCallback`/`ReportCallback` pair stays in place unless a caller (e.g. a benchmark)
asks for `()`.
"""
function build_fluxnet_sim(
    land_domain,
    forcing,
    site,
    toml_dict;
    user_callbacks = nothing,
    diagnostics = nothing,
)
    FT = eltype(land_domain)
    prognostic_land_components = (:canopy, :snow, :soil, :soilco2)
    surface_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)
    soil_domain = land_domain
    (; atmos, radiation) = forcing
    (; site_ID, dt, start_date, stop_date, time_offset, RAI) = site
    (;
        soil_ν,
        soil_K_sat,
        soil_S_s,
        soil_vg_n,
        soil_vg_α,
        θ_r,
        ν_ss_quartz,
        ν_ss_om,
        ν_ss_gravel,
        z_0m_soil,
        z_0b_soil,
        soil_ϵ,
        soil_α_PAR,
        soil_α_NIR,
        Ω,
        χl,
        α_PAR_leaf,
        τ_PAR_leaf,
        α_NIR_leaf,
        τ_NIR_leaf,
        ϵ_canopy,
        ac_canopy,
        g1,
        Vcmax25,
        SAI,
        conductivity_model,
        retention_model,
        plant_ν,
        plant_S_s,
        rooting_depth,
        h_canopy,
    ) = site

    soil = Soil.EnergyHydrology{FT}(
        soil_domain,
        forcing,
        toml_dict;
        prognostic_land_components,
        additional_sources = (ClimaLand.RootExtraction{FT}(),),
        albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(;
            PAR_albedo = soil_α_PAR,
            NIR_albedo = soil_α_NIR,
        ),
        runoff = ClimaLand.Soil.Runoff.SurfaceRunoff(),
        retention_parameters = (;
            ν = soil_ν,
            θ_r,
            K_sat = soil_K_sat,
            hydrology_cm = vanGenuchten{FT}(; α = soil_vg_α, n = soil_vg_n),
        ),
        composition_parameters = (; ν_ss_om, ν_ss_quartz, ν_ss_gravel),
        S_s = soil_S_s,
        z_0m = z_0m_soil,
        z_0b = z_0b_soil,
        emissivity = soil_ϵ,
    )

    co2_prognostic_soil = Soil.Biogeochemistry.PrognosticMet(soil.parameters)
    co2_drivers = Soil.Biogeochemistry.SoilDrivers(co2_prognostic_soil, atmos)
    soilco2 = Soil.Biogeochemistry.SoilCO2Model{FT}(
        soil_domain,
        co2_drivers,
        toml_dict,
    )

    radiative_transfer = Canopy.TwoStreamModel{FT}(
        surface_domain,
        toml_dict;
        radiation_parameters = (;
            Ω,
            G_Function = CLMGFunction(χl),
            α_PAR_leaf,
            τ_PAR_leaf,
            α_NIR_leaf,
            τ_NIR_leaf,
        ),
        ϵ_canopy,
    )
    conductance =
        Canopy.MedlynConductanceModel{FT}(surface_domain, toml_dict; g1)
    photosynthesis = FarquharModel{FT}(
        surface_domain,
        toml_dict;
        photosynthesis_parameters = (; fractional_c3 = FT(1), Vcmax25),
    )

    LAI = ClimaLand.Canopy.prescribed_lai_modis(
        land_domain.space.surface,
        start_date,
        stop_date,
    )
    hydraulics = Canopy.PlantHydraulicsModel{FT}(
        surface_domain,
        toml_dict;
        ν = plant_ν,
        S_s = plant_S_s,
        conductivity_model,
        retention_model,
    )
    biomass = Canopy.PrescribedBiomassModel{FT}(;
        LAI,
        SAI,
        RAI,
        rooting_depth,
        height = h_canopy,
    )
    energy = Canopy.BigLeafEnergyModel{FT}(toml_dict; ac_canopy)
    ground = ClimaLand.PrognosticGroundConditions{FT}()
    canopy = Canopy.CanopyModel{FT}(
        surface_domain,
        (; atmos, radiation, ground),
        LAI,
        toml_dict;
        prognostic_land_components,
        radiative_transfer,
        photosynthesis,
        conductance,
        hydraulics,
        energy,
        biomass,
    )

    snow = Snow.SnowModel(
        FT,
        surface_domain,
        forcing,
        toml_dict,
        dt;
        prognostic_land_components,
    )

    land = LandModel{FT}(canopy, snow, soil, soilco2, nothing)
    set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
        site_ID,
        start_date,
        time_offset,
        land,
    )
    callback_kwarg = isnothing(user_callbacks) ? (;) : (; user_callbacks)
    return LandSimulation(
        start_date,
        stop_date,
        dt,
        land;
        set_ic!,
        updateat = dt,
        diagnostics,
        callback_kwarg...,
    )
end

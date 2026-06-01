###################################
#        MODULE FUNCTIONS         #
###################################

"""
    FluxnetSimulations.generic_site_simulation(site_ID; kwargs...)

Build (but do not solve) a single-column ClimaLand simulation at any fluxnet
site, using default parameters at the site coordinates and MODIS LAI
climatology.

If `lat`, `long`, `time_offset`, `atmos_h` are not provided, they are looked up
from the FLUXNET2015 metadata CSV via [`get_site_info`](@ref) — this requires
the `fluxnet2015` artifact, which is HPC-only by default. On machines without
that artifact (e.g. local development with only the `fluxnet_sites` artifact),
pass coordinates explicitly. For the four bundled sites, you can do this by
calling `get_location(FT, Val(replace_hyphen(site_ID)))` and
`get_fluxtower_height(FT, Val(replace_hyphen(site_ID)))`.

Soil and snow components use the high-level `LandModel` defaults, which fetch
spatially-varying parameters at the column's `(long, lat)`. The canopy is a
PModel + PModelConductance + PiecewiseMoistureStress + PrescribedBiomass
combination, mirroring `experiments/calibration/models/snowy_land.jl` so that
the same calibration knobs (e.g. `pmodel_cstar`, `moisture_stress_c`) work for
both global and single-site runs. Pass `canopy = ...` to override.

# Returns

A NamedTuple `(; simulation, diags, land_domain, start_date, stop_date)`. The
caller is responsible for `solve!(simulation)`.

# Example

    using Dates
    import ClimaLand
    import ClimaLand.FluxnetSimulations as FluxnetSimulations
    using ClimaLand.Simulations: solve!

    # On HPC with fluxnet2015 staged, coordinates auto-resolve from metadata:
    result = FluxnetSimulations.generic_site_simulation("DE-Tha"; duration = Day(7))
    solve!(result.simulation)

    # Locally, pass coordinates from the per-site dispatcher:
    site_val = FluxnetSimulations.replace_hyphen("US-MOz")
    (; lat, long, time_offset) = FluxnetSimulations.get_location(Float64, Val(site_val))
    (; atmos_h) = FluxnetSimulations.get_fluxtower_height(Float64, Val(site_val))
    result = FluxnetSimulations.generic_site_simulation(
        "US-MOz"; lat, long, time_offset, atmos_h, duration = Day(7),
    )
    solve!(result.simulation)
"""
function FluxnetSimulations.generic_site_simulation(
    site_ID::AbstractString;
    FT::Type{<:AbstractFloat} = Float64,
    kwargs...,
)
    # Dispatch through a type-parameterized inner function so `FT` is
    # specialized at the type level — kwargs alone are not specialized in
    # Julia, which produces "Unreachable reached" code-gen failures inside
    # parametric constructors like `PModel{FT}(...)`.
    return _generic_site_simulation(FT, site_ID; kwargs...)
end

function _generic_site_simulation(
    ::Type{FT},
    site_ID::AbstractString;
    lat = nothing,
    long = nothing,
    time_offset = nothing,
    atmos_h = nothing,
    duration::Union{Nothing, Dates.Period} = nothing,
    start_offset::Dates.Period = Dates.Second(0),
    dt = 450.0,
    diagnostic_period = :halfhourly,
    output_vars = [
        "gpp",
        "lhf",
        "shf",
        "lwu",
        "swu",
        "swc",
        "tsoil",
        "swe",
        "et",
    ],
    output_writer = nothing,
    outdir = nothing,
    domain_kwargs = (;
        dz_bottom = 1.5,
        dz_top = 0.1,
        nelements = 20,
        zmin = -10.0,
        zmax = 0.0,
    ),
    toml_dict = LP.create_toml_dict(FT),
    prognostic_land_components = (:canopy, :snow, :soil, :soilco2),
    canopy = nothing,
    soil = nothing,
    snow = nothing,
    space_type = "single",
) where {FT <: AbstractFloat}
    # 1. Resolve coordinates from FLUXNET2015 metadata if any are missing
    if any(isnothing, (lat, long, time_offset, atmos_h))
        info = FluxnetSimulations.get_site_info(site_ID)
        lat = isnothing(lat) ? FT(info.lat) : FT(lat)
        long = isnothing(long) ? FT(info.long) : FT(long)
        time_offset =
            isnothing(time_offset) ? Int(info.time_offset) : time_offset
        atmos_h = if isnothing(atmos_h)
            isempty(info.atmospheric_sensor_height) && error(
                "No atmospheric sensor height in metadata for $site_ID; pass `atmos_h` kwarg.",
            )
            FT(first(info.atmospheric_sensor_height))
        else
            FT(atmos_h)
        end
    else
        lat, long, atmos_h = FT(lat), FT(long), FT(atmos_h)
    end

    # 2. Column domain at (long, lat) — enables spatially-varying defaults
    if space_type == "single"
    land_domain = ClimaLand.Domains.Column(;
        zlim = (FT(domain_kwargs.zmin), FT(domain_kwargs.zmax)),
        nelements = domain_kwargs.nelements,
        dz_tuple = (FT(domain_kwargs.dz_bottom), FT(domain_kwargs.dz_top)),
        longlat = (long, lat),
    )
    else
    land_domain = ClimaLand.Domains.ColumnEnsemble(;
        zlim = (FT(domain_kwargs.zmin), FT(domain_kwargs.zmax)),
        nelements = domain_kwargs.nelements,
        longlat = (long, lat),
        dz_tuple = (FT(domain_kwargs.dz_bottom), FT(domain_kwargs.dz_top)),
    )
    end

    surface_domain = ClimaLand.Domains.obtain_surface_domain(land_domain)
    surface_space = land_domain.space.surface

    # 3. Date range from the site CSV. Pass `required_columns` so the start
    # is advanced past any leading rows where a forcing variable is missing —
    # otherwise the resulting TVIs in step 4 may not cover simulation t=0.
    (start_date, stop_date) = FluxnetSimulations.get_data_dates(
        site_ID,
        time_offset;
        duration,
        start_offset,
        required_columns = FLUXNET_FORCING_COLUMNS,
    )

    # 4. Atmospheric & radiative forcing
    (; atmos, radiation) = FluxnetSimulations.prescribed_forcing_fluxnet(
        site_ID,
        lat,
        long,
        time_offset,
        atmos_h,
        start_date,
        toml_dict,
        FT,
    )
    forcing = (; atmos, radiation)

    # 5. LAI from MODIS climatology (year-agnostic, periodic)
    LAI = ClimaLand.Canopy.prescribed_climatological_lai_modis(surface_space)
    # Main.@infiltrate

    # 6. PModel canopy by default; user can pass `canopy = ...` to override
    if isnothing(canopy)
        photosynthesis = ClimaLand.Canopy.PModel{FT}(land_domain, toml_dict)
        conductance = ClimaLand.Canopy.PModelConductance{FT}(toml_dict)
        soil_moisture_stress =
            ClimaLand.Canopy.PiecewiseMoistureStressModel{FT}(
                land_domain,
                toml_dict,
            )
        biomass = ClimaLand.Canopy.PrescribedBiomassModel{FT}(
            land_domain,
            LAI,
            toml_dict,
        )
        canopy_forcing = (;
            atmos,
            radiation,
            ground = ClimaLand.PrognosticGroundConditions{FT}(),
        )
        canopy = ClimaLand.Canopy.CanopyModel{FT}(
            surface_domain,
            canopy_forcing,
            LAI,
            toml_dict;
            prognostic_land_components,
            photosynthesis,
            conductance,
            soil_moisture_stress,
            biomass,
        )
    end

    # 7. LandModel — uses spatial defaults at longlat for soil, snow, soilco2
    component_kwargs = (; canopy)
    !isnothing(soil) && (component_kwargs = merge(component_kwargs, (; soil)))
    !isnothing(snow) && (component_kwargs = merge(component_kwargs, (; snow)))
    land = ClimaLand.LandModel{FT}(
        forcing,
        LAI,
        toml_dict,
        land_domain,
        dt;
        prognostic_land_components,
        component_kwargs...,
    )

    # 8. Initial conditions from observations at start_date
    set_ic! = FluxnetSimulations.make_set_fluxnet_initial_conditions(
        site_ID,
        start_date,
        time_offset,
        land,
    )

    # 9. Diagnostics
    outdir_resolved = isnothing(outdir) ? mktempdir() : outdir
    output_writer_resolved = if isnothing(output_writer)
        ClimaDiagnostics.Writers.DictWriter()
    else
        output_writer
    end
    diags = ClimaLand.default_diagnostics(
        land,
        start_date,
        outdir_resolved;
        output_writer = output_writer_resolved,
        output_vars,
        reduction_period = diagnostic_period,
    )

    # 10. Build (don't solve — caller decides when)
    simulation = ClimaLand.Simulations.LandSimulation(
        start_date,
        stop_date,
        dt,
        land;
        outdir = outdir_resolved,
        set_ic!,
        updateat = dt,
        diagnostics = diags,
    )

    return (; simulation, diags, land_domain, start_date, stop_date)
end

###################################
#            UTILITIES            #
###################################

"""
    replace_hyphen(old_site_ID::String)

Replaces all instances of hyphens in a given site ID string with underscores
and returns a Symbol of the reformatted site ID to be used as a Val{} type.

For example, an input string "US-MOz" would be output as "US_MOz".
"""
function FluxnetSimulations.replace_hyphen(old_site_ID::String)
    new_site_ID = replace(old_site_ID, "-" => "_")

    return Symbol(new_site_ID)
end

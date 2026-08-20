# Base.show and Base.summary methods for ClimaLand types, collected here
# (rather than beside each type definition) since they are purely
# presentational and span many submodules. Included last in ClimaLand.jl so
# that every type/submodule referenced below has already been defined.

import ClimaComms

### Helpers

# The `space` and (when present) `fields` fields hold the ClimaCore spaces and
# coordinate fields, whose default show output is enormous; every other field
# is a small scalar/tuple that fully describes the domain's configuration.
# This generic fallback is a safety net for any future AbstractDomain subtype
# that doesn't get its own more informative method below.
_domain_show_fieldnames(domain::AbstractDomain) =
    filter(f -> f ∉ (:space, :fields), fieldnames(typeof(domain)))

_domain_device_str(domain::AbstractDomain) =
    nameof(typeof(ClimaComms.device(domain)))

# prognostic_vars can be a flat tuple of symbols (e.g. soil) or a NamedTuple
# grouping symbols by sub-component with many empty entries (e.g. canopy);
# flatten to just the active variable names either way.
_flatten_symbols(x::Symbol) = (x,)
_flatten_symbols(x::Union{Tuple, NamedTuple}) =
    isempty(x) ? () : Tuple(Iterators.flatten(map(_flatten_symbols, values(x))))

# Many parameters (porosity, hydraulic conductivity, leaf reflectance, ...)
# may be either a scalar or a non-scalar value (a ClimaCore Field varying in
# space, or a TimeVaryingInput varying in time); non-scalar values have their
# own default show that's far too long to embed inline in a parameter
# struct's show method, so print a value for scalars and just the type name
# otherwise.
_scalar_or_field_str(x::AbstractFloat) = string(x)
_scalar_or_field_str(x) = string(nameof(typeof(x)))

# A closure/parameterization choice (e.g. a soil hydrology closure) may be
# stored as a scalar struct, or as a Field of that struct when it varies in
# space; either way, name the closure itself rather than saying "Field".
_kind_str(x) = string(nameof(typeof(x)))
_kind_str(x::ClimaCore.Fields.Field) = string(nameof(eltype(x)))

# The precipitation/temperature/wind/humidity/pressure/co2 fields hold
# TimeVaryingInputs backed by NetCDF readers, whose default show output can
# run into the millions of characters; print only the kind of driver
# (e.g. data-backed vs. analytic) instead of the driver's own contents.
_driver_kind(x) = nameof(typeof(x))

### Domains

function Base.show(io::IO, ::MIME"text/plain", domain::AbstractDomain)
    if get(io, :compact, false)
        show(io, domain)
    else
        println(io, nameof(typeof(domain)), "{", eltype(domain), "}")
        for f in _domain_show_fieldnames(domain)
            println(io, "  ", f, ": ", getfield(domain, f))
        end
        println(io, "  device: ", _domain_device_str(domain))
    end
end

function Base.show(io::IO, domain::AbstractDomain)
    print(io, nameof(typeof(domain)), "{", eltype(domain), "}(")
    fs = _domain_show_fieldnames(domain)
    print(io, join(("$f=$(getfield(domain, f))" for f in fs), ", "))
    print(io, ")")
end

Base.summary(io::IO, domain::AbstractDomain) = show(io, domain)

function Base.show(io::IO, ::MIME"text/plain", domain::Point)
    if get(io, :compact, false)
        show(io, domain)
    else
        println(io, "Point{", eltype(domain), "}")
        println(io, "  z_sfc: ", domain.z_sfc, " m")
        println(io, "  device: ", _domain_device_str(domain))
    end
end

Base.show(io::IO, domain::Point) =
    print(io, "Point{", eltype(domain), "}(z_sfc=", domain.z_sfc, " m)")

Base.summary(io::IO, domain::Point) = show(io, domain)

function Base.show(io::IO, ::MIME"text/plain", domain::AbstractColumnDomain)
    if get(io, :compact, false)
        show(io, domain)
    else
        println(io, nameof(typeof(domain)), "{", eltype(domain), "}")
        println(
            io,
            "  z ∈ [",
            domain.zlim[1],
            ", ",
            domain.zlim[2],
            "] m, ",
            domain.nelements[1],
            " elements",
        )
        println(io, "  boundary names: ", domain.boundary_names)
        if !isnothing(domain.dz_tuple)
            println(
                io,
                "  mesh stretching: (dz_bottom=",
                domain.dz_tuple[1],
                ", dz_top=",
                domain.dz_tuple[2],
                ") m",
            )
        end
        println(io, "  device: ", _domain_device_str(domain))
    end
end

Base.show(io::IO, domain::AbstractColumnDomain) = print(
    io,
    nameof(typeof(domain)),
    "{",
    eltype(domain),
    "}(z ∈ [",
    domain.zlim[1],
    ", ",
    domain.zlim[2],
    "], ",
    domain.nelements[1],
    " elements)",
)

Base.summary(io::IO, domain::AbstractColumnDomain) = show(io, domain)

function Base.show(io::IO, ::MIME"text/plain", domain::Plane)
    if get(io, :compact, false)
        show(io, domain)
    else
        println(io, "Plane{", eltype(domain), "}")
        println(
            io,
            "  x ∈ [",
            domain.xlim[1],
            ", ",
            domain.xlim[2],
            "], y ∈ [",
            domain.ylim[1],
            ", ",
            domain.ylim[2],
            "]",
        )
        println(
            io,
            "  nelements: ",
            domain.nelements,
            ", npolynomial: ",
            domain.npolynomial,
            ", periodic: ",
            domain.periodic,
        )
        if !isnothing(domain.longlat)
            println(io, "  longlat center: ", domain.longlat)
        end
        println(io, "  device: ", _domain_device_str(domain))
    end
end

Base.show(io::IO, domain::Plane) = print(
    io,
    "Plane{",
    eltype(domain),
    "}(",
    domain.nelements[1],
    "×",
    domain.nelements[2],
    " elements)",
)

Base.summary(io::IO, domain::Plane) = show(io, domain)

function Base.show(io::IO, ::MIME"text/plain", domain::HybridBox)
    if get(io, :compact, false)
        show(io, domain)
    else
        println(io, "HybridBox{", eltype(domain), "}")
        println(
            io,
            "  x ∈ [",
            domain.xlim[1],
            ", ",
            domain.xlim[2],
            "], y ∈ [",
            domain.ylim[1],
            ", ",
            domain.ylim[2],
            "], z ∈ [",
            domain.zlim[1],
            ", ",
            domain.zlim[2],
            "]",
        )
        println(
            io,
            "  nelements: ",
            domain.nelements,
            ", npolynomial: ",
            domain.npolynomial,
            ", periodic: ",
            domain.periodic,
        )
        if !isnothing(domain.longlat)
            println(io, "  longlat center: ", domain.longlat)
        end
        if !isnothing(domain.dz_tuple)
            println(
                io,
                "  mesh stretching: (dz_bottom=",
                domain.dz_tuple[1],
                ", dz_top=",
                domain.dz_tuple[2],
                ")",
            )
        end
        println(io, "  device: ", _domain_device_str(domain))
    end
end

Base.show(io::IO, domain::HybridBox) = print(
    io,
    "HybridBox{",
    eltype(domain),
    "}(",
    domain.nelements[1],
    "×",
    domain.nelements[2],
    "×",
    domain.nelements[3],
    " elements)",
)

Base.summary(io::IO, domain::HybridBox) = show(io, domain)

function Base.show(io::IO, ::MIME"text/plain", domain::SphericalShell)
    if get(io, :compact, false)
        show(io, domain)
    else
        println(io, "SphericalShell{", eltype(domain), "}")
        println(
            io,
            "  radius: ",
            domain.radius,
            " m, depth: ",
            domain.depth,
            " m",
        )
        println(
            io,
            "  nelements: ",
            domain.nelements,
            " (horizontal, vertical), npolynomial: ",
            domain.npolynomial,
        )
        if !isnothing(domain.dz_tuple)
            println(
                io,
                "  mesh stretching: (dz_bottom=",
                domain.dz_tuple[1],
                ", dz_top=",
                domain.dz_tuple[2],
                ") m",
            )
        end
        println(io, "  device: ", _domain_device_str(domain))
    end
end

Base.show(io::IO, domain::SphericalShell) = print(
    io,
    "SphericalShell{",
    eltype(domain),
    "}(radius=",
    domain.radius,
    " m, depth=",
    domain.depth,
    " m, ",
    domain.nelements[1],
    "×",
    domain.nelements[2],
    " elements)",
)

Base.summary(io::IO, domain::SphericalShell) = show(io, domain)

function Base.show(io::IO, ::MIME"text/plain", domain::SphericalSurface)
    if get(io, :compact, false)
        show(io, domain)
    else
        println(io, "SphericalSurface{", eltype(domain), "}")
        println(io, "  radius: ", domain.radius, " m")
        println(
            io,
            "  nelements: ",
            domain.nelements,
            ", npolynomial: ",
            domain.npolynomial,
        )
        println(io, "  device: ", _domain_device_str(domain))
    end
end

Base.show(io::IO, domain::SphericalSurface) = print(
    io,
    "SphericalSurface{",
    eltype(domain),
    "}(radius=",
    domain.radius,
    " m, ",
    domain.nelements,
    " elements)",
)

Base.summary(io::IO, domain::SphericalSurface) = show(io, domain)

### Parameters

function Base.show(io::IO, ::MIME"text/plain", ps::LP.LandParameters)
    if get(io, :compact, false)
        show(io, ps)
    else
        n_scalar = fieldcount(typeof(ps)) - 3
        println(
            io,
            "LandParameters{",
            eltype(ps),
            "} (",
            n_scalar,
            " physical constants)",
        )
        println(io, "  thermo_params: ", nameof(typeof(ps.thermo_params)))
        println(io, "  surf_flux_params: ", nameof(typeof(ps.surf_flux_params)))
        println(io, "  insol_params: ", nameof(typeof(ps.insol_params)))
    end
end

Base.show(io::IO, ps::LP.LandParameters) =
    print(io, "LandParameters{", eltype(ps), "}")

Base.summary(io::IO, ps::LP.LandParameters) = show(io, ps)

### Boundary conditions

function Base.show(io::IO, ::MIME"text/plain", bc::AbstractBC)
    if get(io, :compact, false)
        show(io, bc)
    else
        println(io, nameof(typeof(bc)))
        for f in fieldnames(typeof(bc))
            v = getfield(bc, f)
            if v isa Function
                println(io, "  ", f, ": <prescribed function>")
            else
                println(io, "  ", f, ": ", sprint(show, v))
            end
        end
    end
end

Base.show(io::IO, bc::AbstractBC) = print(io, nameof(typeof(bc)))

Base.summary(io::IO, bc::AbstractBC) = show(io, bc)

### Models

function Base.show(io::IO, ::MIME"text/plain", model::AbstractModel)
    if get(io, :compact, false)
        show(io, model)
    else
        println(
            io,
            nameof(typeof(model)),
            "{",
            typeof(model).parameters[1],
            "}",
        )
        if hasfield(typeof(model), :domain)
            println(io, "  domain: ", sprint(show, model.domain))
        end
        pv = _flatten_symbols(prognostic_vars(model))
        if !isempty(pv)
            println(io, "  prognostic variables: ", pv)
        end
    end
end

Base.show(io::IO, model::AbstractModel) =
    print(io, nameof(typeof(model)), "{", typeof(model).parameters[1], "}")

Base.summary(io::IO, model::AbstractModel) = show(io, model)

### Drivers

function Base.show(io::IO, ::MIME"text/plain", atmos::PrescribedAtmosphere)
    if get(io, :compact, false)
        show(io, atmos)
    else
        println(
            io,
            nameof(typeof(atmos)),
            "{",
            typeof(atmos).parameters[1],
            "}",
        )
        println(io, "  start_date: ", atmos.start_date)
        println(
            io,
            "  reference height: ",
            atmos.h,
            " m, gustiness: ",
            atmos.gustiness,
            " m/s",
        )
        driver_fields = (:liquid_precip, :snow_precip, :T, :u, :q, :P)
        kinds = unique(_driver_kind(getfield(atmos, f)) for f in driver_fields)
        if length(kinds) == 1
            println(
                io,
                "  drivers: ",
                kinds[1],
                " (liquid/snow precip, T, u, q, P)",
            )
        else
            println(io, "  drivers: ", join(kinds, ", "))
        end
        println(io, "  CO2: ", _driver_kind(atmos.c_co2))
    end
end

function Base.show(io::IO, atmos::PrescribedAtmosphere)
    print(
        io,
        nameof(typeof(atmos)),
        "{",
        typeof(atmos).parameters[1],
        "}(start_date=",
        atmos.start_date,
        ")",
    )
end

Base.summary(io::IO, atmos::PrescribedAtmosphere) = show(io, atmos)

function Base.show(
    io::IO,
    ::MIME"text/plain",
    radiation::PrescribedRadiativeFluxes,
)
    if get(io, :compact, false)
        show(io, radiation)
    else
        println(
            io,
            nameof(typeof(radiation)),
            "{",
            typeof(radiation).parameters[1],
            "}",
        )
        println(io, "  start_date: ", radiation.start_date)
        sw_kind = _driver_kind(radiation.SW_d)
        lw_kind = _driver_kind(radiation.LW_d)
        if sw_kind == lw_kind
            println(io, "  drivers (SW_d, LW_d): ", sw_kind)
        else
            println(io, "  drivers: SW_d=", sw_kind, ", LW_d=", lw_kind)
        end
        println(
            io,
            "  cosθs: ",
            isnothing(radiation.cosθs) ? "computed from date and location" :
            "prescribed",
        )
        println(
            io,
            "  frac_diff: ",
            isnothing(radiation.frac_diff) ? "computed empirically" :
            string("prescribed (", _driver_kind(radiation.frac_diff), ")"),
        )
    end
end

function Base.show(io::IO, radiation::PrescribedRadiativeFluxes)
    print(
        io,
        nameof(typeof(radiation)),
        "{",
        typeof(radiation).parameters[1],
        "}(start_date=",
        radiation.start_date,
        ")",
    )
end

Base.summary(io::IO, radiation::PrescribedRadiativeFluxes) = show(io, radiation)

### LandModel

# A LandModel's components hold references to the same atmos/radiation
# drivers (via boundary conditions), so listing just the component names is
# far more useful than the generic AbstractModel field dump.
function Base.show(io::IO, ::MIME"text/plain", land::LandModel)
    if get(io, :compact, false)
        show(io, land)
    else
        println(io, "LandModel{", typeof(land).parameters[1], "}")
        println(io, "  domain: ", sprint(show, land.soil.domain))
        components =
            filter(f -> !isnothing(getfield(land, f)), fieldnames(typeof(land)))
        for (i, f) in enumerate(components)
            branch = i == length(components) ? "└──" : "├──"
            println(io, branch, " ", f, ": ", nameof(typeof(getfield(land, f))))
        end
    end
end

Base.show(io::IO, land::LandModel) =
    print(io, "LandModel{", typeof(land).parameters[1], "}")

Base.summary(io::IO, land::LandModel) = show(io, land)

### Bucket

function Base.show(
    io::IO,
    ::MIME"text/plain",
    m::Bucket.AbstractBucketAlbedoModel,
)
    if get(io, :compact, false)
        show(io, m)
    else
        println(io, nameof(typeof(m)), "{", typeof(m).parameters[1], "}")
        for f in fieldnames(typeof(m))
            println(io, "  ", f, ": ", _scalar_or_field_str(getfield(m, f)))
        end
    end
end

Base.show(io::IO, m::Bucket.AbstractBucketAlbedoModel) =
    print(io, nameof(typeof(m)), "{", typeof(m).parameters[1], "}")

Base.summary(io::IO, m::Bucket.AbstractBucketAlbedoModel) = show(io, m)

function Base.show(io::IO, ::MIME"text/plain", ps::Bucket.BucketModelParameters)
    if get(io, :compact, false)
        show(io, ps)
    else
        println(io, "BucketModelParameters{", typeof(ps).parameters[1], "}")
        println(
            io,
            "  bucket capacity (W_f): ",
            ps.W_f,
            " m, snow melt timescale (τc): ",
            ps.τc,
            " s",
        )
        println(io, "  albedo: ", sprint(show, ps.albedo))
        println(io, "  earth_param_set: ", sprint(show, ps.earth_param_set))
    end
end

Base.show(io::IO, ps::Bucket.BucketModelParameters) =
    print(io, "BucketModelParameters{", typeof(ps).parameters[1], "}")

Base.summary(io::IO, ps::Bucket.BucketModelParameters) = show(io, ps)

### Snow

function Base.show(io::IO, ::MIME"text/plain", ps::SnowParameters)
    if get(io, :compact, false)
        show(io, ps)
    else
        println(io, "SnowParameters{", typeof(ps).parameters[1], "}")
        println(io, "  timestep (Δt): ", ps.Δt, " s")
        println(
            io,
            "  density: ",
            nameof(typeof(ps.density)),
            ", albedo: ",
            nameof(typeof(ps.α_snow)),
        )
        println(
            io,
            "  thermal conductivity: ",
            nameof(typeof(ps.κ_snow)),
            ", cover fraction: ",
            nameof(typeof(ps.scf)),
        )
        println(io, "  surface temperature: ", nameof(typeof(ps.surf_temp)))
        println(io, "  earth_param_set: ", sprint(show, ps.earth_param_set))
    end
end

Base.show(io::IO, ps::SnowParameters) =
    print(io, "SnowParameters{", typeof(ps).parameters[1], "}")

Base.summary(io::IO, ps::SnowParameters) = show(io, ps)

### Soil

function Base.show(io::IO, ::MIME"text/plain", ps::RichardsParameters)
    if get(io, :compact, false)
        show(io, ps)
    else
        println(io, "RichardsParameters{", eltype(ps.ν), "}")
        println(
            io,
            "  porosity (ν): ",
            _scalar_or_field_str(ps.ν),
            ", K_sat: ",
            _scalar_or_field_str(ps.K_sat),
        )
        println(
            io,
            "  S_s: ",
            _scalar_or_field_str(ps.S_s),
            ", θ_r: ",
            _scalar_or_field_str(ps.θ_r),
        )
        println(io, "  hydrology closure: ", _kind_str(ps.hydrology_cm))
    end
end

Base.show(io::IO, ps::RichardsParameters) =
    print(io, "RichardsParameters{", eltype(ps.ν), "}")

Base.summary(io::IO, ps::RichardsParameters) = show(io, ps)

function Base.show(io::IO, ::MIME"text/plain", ps::EnergyHydrologyParameters)
    if get(io, :compact, false)
        show(io, ps)
    else
        println(io, "EnergyHydrologyParameters{", eltype(ps.ν), "}")
        println(
            io,
            "  porosity (ν): ",
            _scalar_or_field_str(ps.ν),
            ", K_sat: ",
            _scalar_or_field_str(ps.K_sat),
        )
        println(
            io,
            "  S_s: ",
            _scalar_or_field_str(ps.S_s),
            ", θ_r: ",
            _scalar_or_field_str(ps.θ_r),
        )
        println(io, "  hydrology closure: ", _kind_str(ps.hydrology_cm))
        println(io, "  albedo: ", nameof(typeof(ps.albedo)))
        println(io, "  earth_param_set: ", sprint(show, ps.earth_param_set))
    end
end

Base.show(io::IO, ps::EnergyHydrologyParameters) =
    print(io, "EnergyHydrologyParameters{", eltype(ps.ν), "}")

Base.summary(io::IO, ps::EnergyHydrologyParameters) = show(io, ps)

### Runoff

function Base.show(
    io::IO,
    ::MIME"text/plain",
    m::Soil.Runoff.AbstractRunoffModel,
)
    if get(io, :compact, false)
        show(io, m)
    else
        println(io, nameof(typeof(m)))
        for f in fieldnames(typeof(m))
            v = getfield(m, f)
            println(io, "  ", f, ": ", _scalar_or_field_str(v))
        end
    end
end

Base.show(io::IO, m::Soil.Runoff.AbstractRunoffModel) =
    print(io, nameof(typeof(m)))

Base.summary(io::IO, m::Soil.Runoff.AbstractRunoffModel) = show(io, m)

Base.show(io::IO, ::MIME"text/plain", s::Soil.Runoff.TOPMODELSubsurfaceRunoff) =
    show(io, s)

Base.show(io::IO, s::Soil.Runoff.TOPMODELSubsurfaceRunoff) = print(
    io,
    "TOPMODELSubsurfaceRunoff{",
    typeof(s).parameters[1],
    "}(R_sb=",
    s.R_sb,
    ", f_over=",
    s.f_over,
    ", explicit=",
    s.explicit,
    ")",
)

Base.summary(io::IO, s::Soil.Runoff.TOPMODELSubsurfaceRunoff) = show(io, s)

### Vegetation

function Base.show(
    io::IO,
    ::MIME"text/plain",
    m::Canopy.AbstractCanopyComponent,
)
    if get(io, :compact, false)
        show(io, m)
    else
        println(io, nameof(typeof(m)), "{", typeof(m).parameters[1], "}")
        if hasfield(typeof(m), :parameters)
            println(io, "  parameters: ", sprint(show, m.parameters))
        end
    end
end

Base.show(io::IO, m::Canopy.AbstractCanopyComponent) =
    print(io, nameof(typeof(m)), "{", typeof(m).parameters[1], "}")

Base.summary(io::IO, m::Canopy.AbstractCanopyComponent) = show(io, m)

function Base.show(io::IO, ::MIME"text/plain", ps::Canopy.TwoStreamParameters)
    if get(io, :compact, false)
        show(io, ps)
    else
        println(io, "TwoStreamParameters{", eltype(ps), "}")
        println(
            io,
            "  PAR: α_leaf=",
            _scalar_or_field_str(ps.α_PAR_leaf),
            ", τ_leaf=",
            _scalar_or_field_str(ps.τ_PAR_leaf),
        )
        println(
            io,
            "  NIR: α_leaf=",
            _scalar_or_field_str(ps.α_NIR_leaf),
            ", τ_leaf=",
            _scalar_or_field_str(ps.τ_NIR_leaf),
        )
        println(
            io,
            "  clumping index (Ω): ",
            _scalar_or_field_str(ps.Ω),
            ", canopy layers: ",
            ps.n_layers,
        )
        println(io, "  leaf angle distribution: ", _kind_str(ps.G_Function))
    end
end

Base.show(io::IO, ps::Canopy.TwoStreamParameters) =
    print(io, "TwoStreamParameters{", eltype(ps), "}")

Base.summary(io::IO, ps::Canopy.TwoStreamParameters) = show(io, ps)

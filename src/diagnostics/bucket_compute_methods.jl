# stored in p

function compute_albedo!(out, Y, p, t, land_model::BucketModel)
    if isnothing(out)
        return copy(p.bucket.α_sfc)
    else
        out .= p.bucket.α_sfc
    end
end

function compute_net_radiation!(out, Y, p, t, land_model::BucketModel)
    if isnothing(out)
        return copy(p.bucket.R_n)
    else
        out .= p.bucket.R_n
    end
end

function compute_surface_temperature!(out, Y, p, t, land_model::BucketModel)
    if isnothing(out)
        return copy(p.bucket.T_sfc)
    else
        out .= p.bucket.T_sfc
    end
end

function compute_surface_specific_humidity!(
    out,
    Y,
    p,
    t,
    land_model::BucketModel,
)
    if isnothing(out)
        return copy(p.bucket.q_sfc)
    else
        out .= p.bucket.q_sfc
    end
end

function compute_latent_heat_flux!(out, Y, p, t, land_model::BucketModel)
    if isnothing(out)
        return copy(p.bucket.turbulent_fluxes.lhf)
    else
        out .= p.bucket.turbulent_fluxes.lhf
    end
end

function compute_aerodynamic_resistance!(out, Y, p, t, land_model::BucketModel)
    if isnothing(out)
        return copy(p.bucket.turbulent_fluxes.r_ae)
    else
        out .= p.bucket.turbulent_fluxes.r_ae
    end
end

function compute_sensible_heat_flux!(out, Y, p, t, land_model::BucketModel)
    if isnothing(out)
        return copy(p.bucket.turbulent_fluxes.shf)
    else
        out .= p.bucket.turbulent_fluxes.shf
    end
end

function compute_vapor_flux!(out, Y, p, t, land_model::BucketModel)
    if isnothing(out)
        return copy(p.bucket.turbulent_fluxes.vapor_flux)
    else
        out .= p.bucket.turbulent_fluxes.vapor_flux
    end
end

function compute_surface_air_density!(out, Y, p, t, land_model::BucketModel)
    if isnothing(out)
        return copy(p.bucket.ρ_sfc)
    else
        out .= p.bucket.ρ_sfc
    end
end

# stored in Y

function compute_soil_temperature!(out, Y, p, t, land_model::BucketModel)
    if isnothing(out)
        return copy(Y.bucket.T)
    else
        out .= Y.bucket.T
    end
end

function compute_subsurface_water_storage!(
    out,
    Y,
    p,
    t,
    land_model::BucketModel,
)
    if isnothing(out)
        return copy(Y.bucket.W)
    else
        out .= Y.bucket.W
    end
end

function compute_surface_water_content!(out, Y, p, t, land_model::BucketModel)
    if isnothing(out)
        return copy(Y.bucket.Ws)
    else
        out .= Y.bucket.Ws
    end
end

function compute_snow_water_equivalent!(out, Y, p, t, land_model::BucketModel)
    if isnothing(out)
        return copy(Y.bucket.σS)
    else
        out .= Y.bucket.σS
    end
end

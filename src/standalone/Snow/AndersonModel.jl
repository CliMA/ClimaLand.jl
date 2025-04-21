"""
    Anderson1976{FT <: AbstractFloat} <: AbstractDensityModel{FT}
Establishes the density parameterization where snow density
compacts according to the seminal works of Anderson in the 1970s (and specifically the
numerical implementation documented in the Snow17 model).
This uses the default constants defined in the Snow17 documentation, but allows for redefinition of the model constants.
For consistent translation from Snow17 runs and literature, constants are preserved in their original units (the associated `Anderson1976` functions
handle the subsequent unit converstions). As Snow17 evolves unitless relative densities ρ = ρ_snow/ρ_water, some of these constants or rates
are also unitless.
The necessary constants include (preserving constants in their original units)
- c1: rate of (unitless) density increases from overhead load at 0ᵒC (0.026 hr⁻¹ cm⁻¹; cm of water-equivalent load)
- c2: compression constant estimated by Kojima in 1967 (21 cm³/g)
- c3: the "fractional settling" rate (the change in unitless density, per hour, representing destructive metamorphism)
      at 0ᵒC for densities less than the critical density `ρ_d` (default value is 0.005 hr⁻¹)
- c4: the coefficient of temperature-dependent adjustments to the (unitless) fractional settling (0.1 ᵒC⁻¹)
- c5: the multiplicative scaling of the settling rate when no water is present in the snowpack (sources vary on whether this factor should be 0 or 1, default is set to 0)
- cx: the coefficient of unitless-density-dependent adjustments to the (unitless) fractional settling (default value is 23)
- ρ_d: the critical (unitless) density `ρ_d` (default is 0.15)
"""
Base.@kwdef struct Anderson1976{FT} <: AbstractDensityModel{FT}
    c1::FT = FT(0.026)
    c2::FT = FT(21.0)
    c3::FT = FT(0.005)
    c4::FT = FT(0.10)
    c5::FT = FT(0)
    cx::FT = FT(23)
    ρ_d::FT = FT(0.15)
end


#Define the additional prognostic variables needed for using this parameterization:
density_prog_vars(m::Anderson1976) = (:Z,)
density_prog_types(m::Anderson1976{FT}) where {FT} = (FT,)
density_prog_names(m::Anderson1976) = (:surface,)


"""
    compact_density(W_i::FT, sliq::FT, ρ::FT, Δt_t::FT, Tsnow::FT, density::Anderson1976) where {FT}
Returns the new compressed density ρ_new from the input ρ under the compaction model defined
by Snow17 for a given `Anderson1976` parameterization configuration.
"""
function compact_density(
    W_i::FT,
    sliq::FT,
    ρ::FT,
    Δt_t::FT,
    Tsnow::FT,
    density::Anderson1976,
)::FT where {FT}
    c5 = (sliq > 0) ? FT(2) : density.c5
    cx = (ρ > density.ρ_d) ? density.cx : FT(0)

    #Snow17 includes 0.1 * W_i with W_i in millimeters, but we have 100 * W_i as we provide W_i in meters.
    B = 100 * W_i * density.c1 * Δt_t * exp(FT(0.08) * Tsnow - density.c2 * ρ)
    factor = (W_i > 0.001) ? (exp(B) - 1) / B : FT(1)
    A = exp(
        density.c3 *
        c5 *
        Δt_t *
        exp(density.c4 * Tsnow - cx * (ρ - density.ρ_d)),
    )
    factor = factor * A

    ρ_new = minimum([ρ * factor, FT(0.45)])
    ρ_new = (ρ_new < FT(0.05)) ? ρ : ρ_new
    return ρ_new
end


"""
    newsnow_density(air_temp::FT)::FT where {FT}
Estimates the density of newly fallen snow as a function of air temperature (in celsius), in relative units (i.e. as a fraction of the density
of water), according to the Snow17 model.
Used in computing the time derivative of snow depth under the Anderson1976 density parameterization.
"""
function newsnow_density(air_temp::FT)::FT where {FT}
    if air_temp < -15
        return FT(0.05)
    else
        return FT(0.05) + FT(0.0017) * FT((air_temp + FT(15.0))^1.5)
    end
end


"""
    newsnow_temp(air_temp::FT)::FT where {FT}
Estimates the temperature of newly fallen snow as a function of air temperature, according to the Snow17 model.
Used in computing the time derivative of snow depth under the Anderson1976 density parameterization.
"""
function newsnow_temp(air_temp::FT)::FT where {FT}
    return (air_temp > 0) ? FT(0) : air_temp
end


#This is already defined once in NeuralDepthModel.jl and identical for this model type.
#=
"""
    swe_snow_area(S::FT, scf::FT, z::FT)::FT where {FT}
Estimates the SWE over the snow-covered portion of a grid cell, assming Z is provided
as the depth over the snow-covered portion and S is provided as the depth over the whole cell area.
Used in computing the time derivative of snow depth.
"""
function swe_snow_area(S::FT, scf::FT, z::FT)::FT where {FT}
    min_scf = 0.05
    return min(z, S/max(scf, min_scf))
end
=#


"""
    dzdt(density::Anderson1976{FT}, model::SnowModel{FT}, Y, p, t) where {FT}
Returns the change in snow depth (rate) given the current model state and the `Anderson1976`
density paramterization, which estimates contributions from new snowfall, as well as compaction
effects from the existing and newly fallen snow.
"""
function get_dzdt(
    density::Anderson1976{FT},
    model::SnowModel{FT},
    Y,
    p,
) where {FT}
    Δt_t = model.parameters.Δt / 3600 #in hours (when Δt in SnowParameters is ::FT instead of ::Period)

    #Contribution from new snowfall: (uses parameterization of newsnow_temp, newsnow_density)
    air_temp = p.drivers.T .- model.parameters.earth_param_set.T_freeze
    snowfall = abs.(p.drivers.P_snow) .* model.parameters.Δt
    T_newsnow = newsnow_temp.(air_temp)
    ρ_precip = newsnow_density.(air_temp) #unitless
    ρ_newsnow =
        compact_density.(
            snowfall,
            FT(0),
            ρ_precip,
            Δt_t,
            T_newsnow,
            Ref(density),
        )
    Δz = snowfall ./ ρ_newsnow

    #Estimate resulting change from compaction effects of old snow (and the incoming snow on top):
    S_snowheight =
        swe_snow_area.(Y.snow.S, p.snow.snow_cover_fraction, Y.snow.Z)
    W_ice = (1 .- p.snow.q_l) .* S_snowheight #S_snowheight instead of Y.snow.S
    ρ_ice = W_ice ./ Y.snow.Z #unitless
    sliq = p.snow.q_l .* S_snowheight
    Tsnow = p.snow.T .- model.parameters.earth_param_set.T_freeze
    W_use = W_ice .+ snowfall
    parent(ρ_ice)[parent(W_ice) .<= eps(FT)] .= FT(0.1) #any nonzero value for no-snowpack to avoid NaN and return 0.0
    ρ_est = compact_density.(W_use, sliq, ρ_ice, Δt_t, Tsnow, Ref(density))
    dz_compaction = (W_ice ./ ρ_est) .- Y.snow.Z

    return (dz_compaction .+ Δz) ./ model.parameters.Δt
end


# Do we want to change the type assertion to ::Union{Anderson1976, NeuralDepthModel} in NeuralDepthModel
# to avoid redundancy since the two are identical and use the same functionality, or is better to keep
# them separate and visible with their own code file? I can always just comment this one out like I did
# above/below with clip_dZdt() and swe_snow_area. 
"""
    update_density_and_depth!(ρ_snow, z_snow,density::Anderson1976, Y, p, params::SnowParameters,)

Updates the snow density and depth in place given the current model state. Extends the update_density_and_depth!
function for the Anderson1976 type.
"""
function update_density_and_depth!(
    ρ_snow,
    z_snow,
    density::Anderson1976,
    Y,
    p,
    params::SnowParameters,
)
    @. z_snow = Y.snow.Z
    @. ρ_snow = snow_bulk_density(Y.snow.S, z_snow, params)
end

#This method already written once in NeuralDepthModel.jl and applies here too:
#=
"""
    clip_dZdt(S::FT, Z::FT, dSdt::FT, dZdt::FT, Δt::FT)::FT

A helper function which clips the tendency of Z such that
its behavior is consistent with that of S: if all snow melts
within a timestep, we clip the tendency of S so that it does
not become negative, and here we also clip the tendency of Z
so that depth does not become negative. Additionally, if the
tendencies of Z and S are such that we would encounter Z < S
(rho_snow > rho_liq, for Z and S over the snow-area, not
the grid cell area), we also clip the tendency.
"""
function clip_dZdt(
    S::FT,
    Z::FT,
    dSdt::FT,
    dZdt::FT,
    scf::FT,
    Δt::FT,
)::FT where {FT}
    min_scf = 0.05
    #Case if S is set to zero:
    if (S + dSdt * Δt) <= eps(FT)
        return -Z / Δt
        #Case if Z would have been set to Z < S:
    elseif (Z + dZdt * Δt) < ((S + dSdt * Δt) / max(scf, eps(FT)))
        #more stable form for very small Z, S
        return ((dSdt * Δt + S) / max(scf, min_scf) - Z) / Δt
    else
        return dZdt
    end
end
=#

"""
    update_density_prog!(density::Anderson1976, model::SnowModel, Y, p)

Updates all prognostic variables associated with density/depth given the current model state and the `Anderson1976`
density paramterization.
"""
function update_density_prog!(density::Anderson1976, model::SnowModel, dY, Y, p)
    dY.snow.Z .= get_dzdt(density, model, Y, p)

    # Now we clip the tendency so that Z stays within approximately physical bounds.
    @. dY.snow.Z = clip_dZdt(
        Y.snow.S,
        Y.snow.Z,
        dY.snow.S, #assumes dY.snow.S is updated (and clipped) before dY.snow.Z
        dY.snow.Z,
        p.snow.snow_cover_fraction,
        model.parameters.Δt,
    )
end

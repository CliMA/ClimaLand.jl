import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.Regridders: InterpolationsRegridder
import Interpolations

export make_set_initial_state_from_file

regridder_type = :InterpolationsRegridder
extrapolation_bc =
    (Interpolations.Periodic(), Interpolations.Flat(), Interpolations.Flat())
"""
    set_soil_initial_conditions!(Y, ν, θ_r, subsurface_space, soil_ic_path)

Sets the soil initial conditions, stored in `Y.soil.ϑ_l`, `Y.soil.θ_i`,
`Y.soil.ρe_int`, and defined on the `subsurface_space`, using the values
in the net cdf file stored at `soil_ic_path`.

Since the values in the netcdf file have been interpolated once (to save the output),
and interpolated again (onto the model grid), there is no guarantee that the liquid water
content is above the residual. Although our model can simulate oversaturated soils, there
is no guarantee that the initial conditions read in will be stable. Because of this,
we enforce the constraint of ϑ_l > θ_r and θ_i + ϑ_l < ν.
"""
function set_soil_initial_conditions!(
    Y,
    ν,
    θ_r,
    subsurface_space,
    soil_ic_path,
    soil,
    T_bounds;
    regridder_type = regridder_type,
    extrapolation_bc = extrapolation_bc,
)
    Y.soil.ϑ_l .= SpaceVaryingInput(
        soil_ic_path,
        "swc",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )
    Y.soil.θ_i .= SpaceVaryingInput(
        soil_ic_path,
        "si",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )

    Y.soil.ϑ_l .= enforce_residual_constraint.(Y.soil.ϑ_l, θ_r)
    Y.soil.ϑ_l .= enforce_porosity_constraint.(Y.soil.ϑ_l, ν)
    Y.soil.θ_i .=
        enforce_residual_constraint.(Y.soil.θ_i, eltype(Y.soil.θ_i)(0))
    Y.soil.θ_i .= enforce_porosity_constraint.(Y.soil.ϑ_l, Y.soil.θ_i, ν)
    ρc_s =
        ClimaLand.Soil.volumetric_heat_capacity.(
            Y.soil.ϑ_l,
            Y.soil.θ_i,
            soil.parameters.ρc_ds,
            soil.parameters.earth_param_set,
        )
    Y.soil.ρe_int .= SpaceVaryingInput(
        soil_ic_path,
        "sie",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )
    T =
        ClimaLand.Soil.temperature_from_ρe_int.(
            Y.soil.ρe_int,
            Y.soil.θ_i,
            ρc_s,
            soil.parameters.earth_param_set,
        )
    T .= clip_to_bounds.(T, T_bounds[1], T_bounds[2])
    Y.soil.ρe_int .=
        ClimaLand.Soil.volumetric_internal_energy.(
            Y.soil.θ_i,
            ρc_s,
            T,
            soil.parameters.earth_param_set,
        )
    return nothing
end

"""
    clip_to_bounds(
    T::FT,
    lb::FT,
    ub::FT,
) where {FT <: AbstractFloat}

Clips `T` to be in the bounds [lb, ub].

In our initial condition function, this is used to ensure that the temperature of the soil
is within a physical range.
"""
function clip_to_bounds(T::FT, lb::FT, ub::FT) where {FT <: AbstractFloat}
    if T > ub
        return ub
    elseif T < lb
        return lb
    else
        return T
    end
end

"""
     enforce_residual_constraint(ϑ_l::FT ,θ_r::FT)

Enforces the constraint that ϑ_l > θ_r by returning 1.05 θ_r
if ϑ_l < θ_r, and ϑ_l otherwise.
"""
function enforce_residual_constraint(ϑ_l::FT, θ_r::FT) where {FT}
    if ϑ_l < θ_r
        return θ_r * FT(1.05)
    else
        return ϑ_l
    end
end


"""
    enforce_porosity_constraint(ϑ_l::FT, ν::FT)

Enforces the constraint that ϑ_l <= ν.
"""
function enforce_porosity_constraint(ϑ_l::FT, ν::FT) where {FT}
    if ϑ_l > ν # if we exceed porosity
        return FT(0.95) * ν # clip water content to 95% of available pore space
    else
        return ϑ_l
    end
end
"""
        enforce_porosity_constraint(ϑ_l::FT, θ_i::FT, ν::FT, θ_r::FT)
    
Enforces the constraint that ϑ_l + θ_i <= ν, by clipping the ice content to be
99% (ν-ϑ_l), or leaving it unchanged if the constraint is already satisfied.
"""
function enforce_porosity_constraint(ϑ_l::FT, θ_i::FT, ν::FT) where {FT}
    if ϑ_l + θ_i > ν # if we exceed porosity
        return FT(0.95) * (ν - ϑ_l) # clip ice content to 95% of available pore space
    else
        return θ_i
    end
end


"""
    set_snow_initial_conditions!(Y, p, surface_space, snow_ic_path, params)

Sets the snow initial conditions, stored in `Y.snow.S`, `Y.snow.S_l`,
`Y.snow.U`, and defined on the `surface_space`, using the values
in the net cdf file stored at `snow_ic_path`.

We assume that S_l = 0 and that the temperature of the snow has been updated
in p.snow.T.
"""
function set_snow_initial_conditions!(
    Y,
    p,
    surface_space,
    snow_ic_path,
    params;
    regridder_type = regridder_type,
    extrapolation_bc = extrapolation_bc,
)
    Y.snow.S .=
        SpaceVaryingInput(
            snow_ic_path,
            "swe",
            surface_space;
            regridder_type,
            regridder_kwargs = (; extrapolation_bc,),
        ) .* 0
    Y.snow.S_l .= 0
    p.snow.T .= enforce_snow_temperature_constraint.(Y.snow.S, p.snow.T)
    Y.snow.U .=
        ClimaLand.Snow.energy_from_T_and_swe.(Y.snow.S, p.snow.T, params)
    return nothing
end
"""
    enforce_snow_temperature_constraint(S::FT, T::FT)

Enforces the constraint that T < 273.15 if S > 0.
"""
function enforce_snow_temperature_constraint(S::FT, T::FT) where {FT}
    if S > sqrt(eps(FT)) # if snow is on the ground
        return min(T, FT(273))
    else
        return T
    end
end

"""
    make_set_initial_state_from_file(ic_path, land::LandModel{FT}) where {FT}

Returns a function which takes (Y,p,t0,land) as arguments, and updates
the state Y in place with initial conditions from `ic_path`, a netCDF file.
Fields in the cache `p` are used as pre-allocated memory and are updated as
well, but this does not mean that the cache state is consitent with Y and t entirely.

Currently only tested and used for global simulations, but the same returned
function should work for column simulations.

The returned function is a closure for `ic_path`. It could also be for `land`, as
many other ClimaLand functions are, but we wish to preserve the argument `land`
in `set_ic!` for users who wish to define their own initial condition function,
which may require parameters, etc, stored in `land`.
"""
function make_set_initial_state_from_file(
    ic_path,
    land::LandModel{FT},
) where {FT}
    function set_ic!(Y, p, t0, land)
        atmos = land.soil.boundary_conditions.top.atmos
        evaluate!(p.snow.T, atmos.T, t0)
        set_snow_initial_conditions!(
            Y,
            p,
            land.snow.domain.space.surface,
            ic_path,
            land.snow.parameters,
        )
        Y.soilco2.C .= FT(0.000412) # set to atmospheric co2, mol co2 per mol air
        Y.canopy.hydraulics.ϑ_l.:1 .= land.canopy.hydraulics.parameters.ν
        evaluate!(Y.canopy.energy.T, atmos.T, t0)
        T_bounds = extrema(Y.canopy.energy.T)

        set_soil_initial_conditions!(
            Y,
            land.soil.parameters.ν,
            land.soil.parameters.θ_r,
            land.soil.domain.space.subsurface,
            ic_path,
            land.soil,
            T_bounds,
        )
    end
    return set_ic!
end

"""
    make_set_initial_state_from_file(ic_path, land::SoilCanopyModel{FT}) where {FT}

Returns a function which takes (Y,p,t0,land) as arguments, and updates
the state Y in place with initial conditions from `ic_path`, a netCDF file.
Fields in the cache `p` are used as pre-allocated memory and are updated as
well, but this does not mean that the cache state is consitent with Y and t entirely.

Currently only tested and used for global simulations, but the same returned
function should work for column simulations.

The returned function is a closure for `ic_path`. It could also be for `land`, as
many other ClimaLand functions are, but we wish to preserve the argument `land`
in `set_ic!` for users who wish to define their own initial condition function,
which may require parameters, etc, stored in `land`.
"""
function make_set_initial_state_from_file(
    ic_path,
    land::SoilCanopyModel{FT},
) where {FT}
    function set_ic!(Y, p, t0, land)
        atmos = land.soil.boundary_conditions.top.atmos
        evaluate!(p.drivers.T, atmos.T, t0)

        Y.soilco2.C .= FT(0.000412) # set to atmospheric co2, mol co2 per mol air
        Y.canopy.hydraulics.ϑ_l.:1 .= land.canopy.hydraulics.parameters.ν
        evaluate!(Y.canopy.energy.T, atmos.T, t0)
        T_bounds = extrema(Y.canopy.energy.T)

        set_soil_initial_conditions!(
            Y,
            land.soil.parameters.ν,
            land.soil.parameters.θ_r,
            land.soil.domain.space.subsurface,
            ic_path,
            land.soil,
            T_bounds,
        )
    end
    return set_ic!
end

"""
    make_set_initial_state_from_file(ic_path, model::ClimaLand.Soil.EnergyHydrology{FT}) where {FT}

Returns a function which takes (Y,p,t0,model) as arguments, and updates
the state Y in place with initial conditions from `ic_path`, a netCDF file.
Fields in the cache `p` are used as pre-allocated memory and are updated as
well, but this does not mean that the cache state is consitent with Y and t entirely.

Currently only tested and used for global simulations, but the same returned
function should work for column simulations.

The returned function is a closure for `ic_path`. It could also be for `model`, as
many other ClimaLand functions are, but we wish to preserve the argument `model`
in `set_ic!` for users who wish to define their own initial condition function,
which may require parameters, etc, stored in `model`.
"""
function make_set_initial_state_from_file(
    ic_path,
    model::ClimaLand.Soil.EnergyHydrology{FT},
) where {FT}
    function set_ic!(Y, p, t0, model)
        atmos = model.boundary_conditions.top.atmos
        evaluate!(p.drivers.T, atmos.T, t0)
        T_bounds = extrema(p.drivers.T)

        set_soil_initial_conditions!(
            Y,
            model.parameters.ν,
            model.parameters.θ_r,
            model.domain.space.subsurface,
            ic_path,
            model,
            T_bounds,
        )
    end
    return set_ic!
end

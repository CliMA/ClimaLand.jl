import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.Regridders: InterpolationsRegridder
import Interpolations
using ClimaCore

export make_set_initial_state_from_file

regridder_type = :InterpolationsRegridder
extrapolation_bc =
    (Interpolations.Periodic(), Interpolations.Flat(), Interpolations.Flat())
interpolation_method = Interpolations.Constant()
"""
    set_soil_initial_conditions!(Y, subsurface_space, soil_ic_path, soil; enforce_constraints = false, T_bounds = nothing)

Sets the soil initial conditions, stored in `Y.soil.ϑ_l`, `Y.soil.θ_i`,
`Y.soil.ρe_int`, and defined on the `subsurface_space`, using the values
in the net cdf file stored at `soil_ic_path`.

Since the values in the netcdf file have been interpolated once (to save the output),
and interpolated again (onto the model grid), there is no guarantee that the liquid water
content is above the residual, or that the total water is below porosity.
 Although our model can simulate oversaturated soils, there
is no guarantee that the initial conditions read in will be stable. Because of this, you
can optionally enforce the constraint of ϑ_l > θ_r and θ_i + ϑ_l < ν
by setting `enforce_constraints = true`. This will also clip the temperature to be
in the range of `T_bounds`.
"""
function set_soil_initial_conditions!(
    Y,
    subsurface_space,
    soil_ic_path,
    soil,
    ;
    regridder_type = regridder_type,
    extrapolation_bc = extrapolation_bc,
    interpolation_method = interpolation_method,
    enforce_constraints = false,
    T_bounds = nothing,
)
    (; ν, θ_r) = soil.parameters

    Y.soil.ϑ_l .= SpaceVaryingInput(
        soil_ic_path,
        "swc",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
    )
    Y.soil.θ_i .= SpaceVaryingInput(
        soil_ic_path,
        "si",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
    )

    Y.soil.ρe_int .= SpaceVaryingInput(
        soil_ic_path,
        "sie",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
    )
    if enforce_constraints
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
        if ~isnothing(T_bounds)
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
        end
    end
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

Enforces the constraint that ϑ_l > θ_r by returning 1.1 θ_r
if ϑ_l < θ_r, and ϑ_l otherwise.
"""
function enforce_residual_constraint(ϑ_l::FT, θ_r::FT) where {FT}
    if ϑ_l < θ_r * FT(1.05)
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
    if ϑ_l > ν * FT(0.95) # if we exceed porosity
        return FT(0.95) * ν # clip water content to 95% of available pore space
    else
        return ϑ_l
    end
end
"""
        enforce_porosity_constraint(ϑ_l::FT, θ_i::FT, ν::FT, θ_r::FT)

Enforces the constraint that ϑ_l + θ_i <= ν, by clipping the ice content to be
95% (ν-ϑ_l), or leaving it unchanged if the constraint is already satisfied.
"""
function enforce_porosity_constraint(ϑ_l::FT, θ_i::FT, ν::FT) where {FT}
    if ϑ_l + θ_i > FT(0.95) * ν # if we exceed porosity
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
    interpolation_method = interpolation_method,
)
    Y.snow.S .= SpaceVaryingInput(
        snow_ic_path,
        "swe",
        surface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
    )
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
        return min(T, FT(273.15))
    else
        return FT(273.16)
    end
end

"""
    make_set_initial_state_from_file(ic_path, land::LandModel{FT}; enforce_constraints=false) where {FT}

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

If `enforce_constraints = true`, we ensure the soil water content is between porosity and the
residual value, and that the temperature is bounded to be within the extrema of the air temperature
at the surface.

It is assumed that in CoupledAtmosphere simulations that `p.drivers.T` has 
been updated already.
"""
function make_set_initial_state_from_file(
    ic_path,
    land::LandModel{FT};
    enforce_constraints = false,
) where {FT}
    function set_ic!(Y, p, t0, land)
        atmos = land.soil.boundary_conditions.top.atmos
        if atmos isa ClimaLand.PrescribedAtmosphere
            evaluate!(p.drivers.T, atmos.T, t0)
        end
        # SoilCO2 IC
        if !isnothing(land.soilco2)
            Y.soilco2.CO2 .= FT(0.000412) # set to atmospheric co2, mol co2 per mol air
            Y.soilco2.O2_f .= FT(0.21)    # atmospheric O2 volumetric fraction
            Y.soilco2.SOC .= FT(5.0)      # default SOC concentration (kg C/m³)
        end
        # Soil IC
        if enforce_constraints
            # get only the values over land
            T_atmos = Array(parent(p.drivers.T))[:]
            T_bounds = extrema(T_atmos[T_atmos .> sqrt(eps(FT))])
        else
            T_bounds = nothing
        end

        set_soil_initial_conditions!(
            Y,
            land.soil.domain.space.subsurface,
            ic_path,
            land.soil;
            T_bounds,
            enforce_constraints,
        )

        # Snow IC
        # Use soil temperature at top to set IC
        soil = land.soil
        ρc_s =
            ClimaLand.Soil.volumetric_heat_capacity.(
                min.(Y.soil.ϑ_l, soil.parameters.ν .- Y.soil.θ_i), # θ_l
                Y.soil.θ_i,
                soil.parameters.ρc_ds,
                soil.parameters.earth_param_set,
            )
        p.soil.T .=
            ClimaLand.Soil.temperature_from_ρe_int.(
                Y.soil.ρe_int,
                Y.soil.θ_i,
                ρc_s,
                soil.parameters.earth_param_set,
            )
        T_sfc = ClimaLand.Domains.top_center_to_surface(p.soil.T)
        p.snow.T .= T_sfc
        set_snow_initial_conditions!(
            Y,
            p,
            land.snow.domain.space.surface,
            ic_path,
            land.snow.parameters,
        )


        # Canopy IC
        # First determine if leaf water potential is in the file. If so, use
        # that to set the IC; otherwise choose steady state with the soil water.
        ds = NCDataset(ic_path, "r")
        variable_names = keys(ds)
        close(ds)
        if land.canopy.hydraulics isa ClimaLand.Canopy.PlantHydraulicsModel
            if "lwp" ∈ variable_names
                ψ_roots = SpaceVaryingInput(
                    ic_path,
                    "lwp",
                    land.canopy.domain.space.surface;
                    regridder_type,
                    regridder_kwargs = (;
                        extrapolation_bc,
                        interpolation_method,
                    ),
                )
                Y.canopy.hydraulics.ϑ_l.:1 .=
                    ClimaLand.Canopy.PlantHydraulics.inverse_water_retention_curve.(
                        land.canopy.hydraulics.parameters.retention_model,
                        ψ_roots,
                        land.canopy.hydraulics.parameters.ν,
                        land.canopy.hydraulics.parameters.S_s,
                    ) .* land.canopy.hydraulics.parameters.ν
            else
                @. p.soil.ψ = ClimaLand.Soil.pressure_head(
                    land.soil.parameters.hydrology_cm,
                    land.soil.parameters.θ_r,
                    Y.soil.ϑ_l,
                    land.soil.parameters.ν - Y.soil.θ_i,
                    land.soil.parameters.S_s,
                )
                ψ_roots =
                    ClimaCore.Fields.zeros(axes(Y.canopy.hydraulics.ϑ_l.:1))
                z = land.soil.domain.fields.z
                tmp = @. ClimaLand.Canopy.root_distribution(
                    z,
                    land.canopy.biomass.rooting_depth,
                ) * p.soil.ψ / land.soil.domain.fields.depth
                ClimaCore.Operators.column_integral_definite!(ψ_roots, tmp)
            end
        end
        if land.canopy.energy isa ClimaLand.Canopy.BigLeafEnergyModel
            Y.canopy.energy.T .= p.drivers.T
        end
    end
    return set_ic!
end

"""
    make_set_initial_state_from_file(ic_path, land::SoilCanopyModel{FT}; enforce_constraints = false) where {FT}

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

If `enforce_constraints = true`, we ensure the soil water content is between porosity and the
residual value, and that the temperature is bounded to be within the extrema of the air temperature
at the surface.

It is assumed that in CoupledAtmosphere simulations that `p.drivers.T` has 
been updated already.
"""
function make_set_initial_state_from_file(
    ic_path,
    land::SoilCanopyModel{FT};
    enforce_constraints = false,
) where {FT}
    function set_ic!(Y, p, t0, land)
        atmos = land.soil.boundary_conditions.top.atmos
        if atmos isa ClimaLand.PrescribedAtmosphere
            evaluate!(p.drivers.T, atmos.T, t0)
        end

        Y.soilco2.CO2 .= FT(0.000412) # set to atmospheric co2, mol co2 per mol air
        Y.soilco2.O2_f .= FT(0.21)    # atmospheric O2 volumetric fraction
        Y.soilco2.SOC .= FT(5.0)      # default SOC concentration (kg C/m³)

        if enforce_constraints
            # get only the values over land
            T_atmos = Array(parent(p.drivers.T))[:]
            T_bounds = extrema(T_atmos[T_atmos .> sqrt(eps(FT))])
        else
            T_bounds = nothing
        end

        set_soil_initial_conditions!(
            Y,
            land.soil.domain.space.subsurface,
            ic_path,
            land.soil;
            T_bounds,
            enforce_constraints,
        )
        # Canopy IC
        # First determine if leaf water potential is in the file. If so, use
        # that to set the IC; otherwise choose steady state with the soil water.
        ds = NCDataset(ic_path, "r")
        variable_names = keys(ds)
        close(ds)
        if land.canopy.hydraulics isa ClimaLand.Canopy.PlantHydraulicsModel
            if "lwp" ∈ variable_names
                ψ_roots = SpaceVaryingInput(
                    ic_path,
                    "lwp",
                    land.canopy.domain.space.surface;
                    regridder_type,
                    regridder_kwargs = (;
                        extrapolation_bc,
                        interpolation_method,
                    ),
                )
                Y.canopy.hydraulics.ϑ_l.:1 .=
                    ClimaLand.Canopy.PlantHydraulics.inverse_water_retention_curve.(
                        land.canopy.hydraulics.parameters.retention_model,
                        ψ_roots,
                        land.canopy.hydraulics.parameters.ν,
                        land.canopy.hydraulics.parameters.S_s,
                    ) .* land.canopy.hydraulics.parameters.ν
            else
                @. p.soil.ψ = ClimaLand.Soil.pressure_head(
                    land.soil.parameters.hydrology_cm,
                    land.soil.parameters.θ_r,
                    Y.soil.ϑ_l,
                    land.soil.parameters.ν - Y.soil.θ_i,
                    land.soil.parameters.S_s,
                )
                ψ_roots =
                    ClimaCore.Fields.zeros(axes(Y.canopy.hydraulics.ϑ_l.:1))
                z = land.soil.domain.fields.z
                tmp = @. ClimaLand.Canopy.root_distribution(
                    z,
                    land.canopy.biomass.rooting_depth,
                ) * p.soil.ψ / land.soil.domain.fields.depth
                ClimaCore.Operators.column_integral_definite!(ψ_roots, tmp)
            end
        end
        if land.canopy.energy isa ClimaLand.Canopy.BigLeafEnergyModel
            Y.canopy.energy.T .= p.drivers.T
        end
    end
    return set_ic!
end

"""
    make_set_initial_state_from_file(ic_path, model::ClimaLand.Soil.EnergyHydrology{FT}; enforce_constraints = false) where {FT}

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

If `enforce_constraints = true`, we ensure the soil water content is between porosity and the
residual value, and that the temperature is bounded to be within the extrema of the air temperature
at the surface.

It is assumed that in CoupledAtmosphere simulations that `p.drivers.T` has 
been updated already.
"""
function make_set_initial_state_from_file(
    ic_path,
    model::ClimaLand.Soil.EnergyHydrology{FT};
    enforce_constraints = false,
) where {FT}
    function set_ic!(Y, p, t0, model)
        atmos = model.boundary_conditions.top.atmos
        if atmos isa ClimaLand.PrescribedAtmosphere
            evaluate!(p.drivers.T, atmos.T, t0)
        end
        if enforce_constraints
            # get only the values over land
            T_atmos = Array(parent(p.drivers.T))[:]
            T_bounds = extrema(T_atmos[T_atmos .> sqrt(eps(FT))])
        else
            T_bounds = nothing
        end
        set_soil_initial_conditions!(
            Y,
            model.domain.space.subsurface,
            ic_path,
            model;
            T_bounds,
            enforce_constraints,
        )
    end
    return set_ic!
end

function set_soil_initial_conditions_from_temperature_and_total_water!(
    Y,
    subsurface_space,
    soil_ic_path,
    soil;
    water_varname = "swvl",
    temp_varname = "stl",
    regridder_type = regridder_type,
    extrapolation_bc = extrapolation_bc,
    interpolation_method = interpolation_method,
)
    total_water_content = SpaceVaryingInput(
        soil_ic_path,
        water_varname,
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
    )
    temperature = SpaceVaryingInput(
        soil_ic_path,
        temp_varname,
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc, interpolation_method),
    )
    # The interpolation of gridded date to field is not mask aware.
    # This results in unphysical values along the coast.
    # To mitigate, we find the bounds of the land values in the data
    # (ocean values assumed set to zero),
    # and clip the interpolated field to lie within those bounds.
    raw_data = NCDataset(soil_ic_path)
    raw_temp_data = raw_data[temp_varname][:]
    FT = eltype(Y.soil.ϑ_l)
    T_bounds = FT.(extrema(raw_temp_data[raw_temp_data .> 0]))
    close(raw_data)
    temperature .= clip_to_bounds.(temperature, T_bounds[1], T_bounds[2])
    (; θ_r, ν, ρc_ds, earth_param_set) = soil.parameters
    _T_freeze = LP.T_freeze(earth_param_set)
    function liquid_soil_water(twc, T, θ_r, ν)
        if T > _T_freeze
            return twc
        else
            return θ_r
        end
    end

    Y.soil.ϑ_l .= liquid_soil_water.(total_water_content, temperature, θ_r, ν)
    Y.soil.θ_i .= total_water_content .- Y.soil.ϑ_l

    Y.soil.ϑ_l .= enforce_residual_constraint.(Y.soil.ϑ_l, θ_r)
    Y.soil.ϑ_l .= enforce_porosity_constraint.(Y.soil.ϑ_l, ν)
    Y.soil.θ_i .=
        enforce_residual_constraint.(Y.soil.θ_i, eltype(Y.soil.θ_i)(0))
    Y.soil.θ_i .= enforce_porosity_constraint.(Y.soil.ϑ_l, Y.soil.θ_i, ν)
    ρc_s =
        ClimaLand.Soil.volumetric_heat_capacity.(
            Y.soil.ϑ_l,
            Y.soil.θ_i,
            ρc_ds,
            earth_param_set,
        )
    Y.soil.ρe_int .=
        ClimaLand.Soil.volumetric_internal_energy.(
            Y.soil.θ_i,
            ρc_s,
            temperature,
            earth_param_set,
        )
    return nothing
end

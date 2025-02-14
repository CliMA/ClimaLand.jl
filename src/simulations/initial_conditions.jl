import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.Regridders: InterpolationsRegridder
import Interpolations
regridder_type = :InterpolationsRegridder
extrapolation_bc =
    (Interpolations.Periodic(), Interpolations.Flat(), Interpolations.Flat())
"""
    set_soil_initial_conditions!(Y, subsurface_space, soil_ic_path)

Sets the soil initial conditions, stored in `Y.soil.ϑ_l`, `Y.soil.θ_i`,
`Y.soil.ρe_int`, and defined on the `subsurface_space`, using the values
in the net cdf file stored at `soil_ic_path`.
"""
function set_soil_initial_conditions!(
    Y,
    subsurface_space,
    soil_ic_path;
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
    Y.soil.ρe_int .= SpaceVaryingInput(
        soil_ic_path,
        "sie",
        subsurface_space;
        regridder_type,
        regridder_kwargs = (; extrapolation_bc,),
    )
    return nothing
end

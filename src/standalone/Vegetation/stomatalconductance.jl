export MedlynConductanceParameters, MedlynConductanceModel, 
         PModelConductanceParameters, PModelConductance

abstract type AbstractStomatalConductanceModel{FT} <:
              AbstractCanopyComponent{FT} end

"""
    MedlynConductanceParameters{FT <: AbstractFloat}

The required parameters for the Medlyn stomatal conductance model.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct MedlynConductanceParameters{
    FT <: AbstractFloat,
    G1 <: Union{FT, ClimaCore.Fields.Field},
}
    "Relative diffusivity of water vapor (unitless)"
    Drel::FT
    "Minimum stomatal conductance mol/m^2/s"
    g0::FT
    "Slope parameter, inversely proportional to the square root of marginal water use efficiency (Pa^{1/2})"
    g1::G1
end

Base.eltype(::MedlynConductanceParameters{FT}) where {FT} = FT

struct MedlynConductanceModel{FT, MCP <: MedlynConductanceParameters{FT}} <:
       AbstractStomatalConductanceModel{FT}
    parameters::MCP
end

function MedlynConductanceModel{FT}(
    parameters::MedlynConductanceParameters{FT},
) where {FT <: AbstractFloat}
    return MedlynConductanceModel{eltype(parameters), typeof(parameters)}(
        parameters,
    )
end

ClimaLand.name(model::AbstractStomatalConductanceModel) = :conductance

ClimaLand.auxiliary_vars(model::MedlynConductanceModel) = (:r_stomata_canopy,)
ClimaLand.auxiliary_types(model::MedlynConductanceModel{FT}) where {FT} = (FT,)
ClimaLand.auxiliary_domain_names(::MedlynConductanceModel) = (:surface,)

"""
    update_canopy_conductance!(p, Y, model::MedlynConductanceModel, canopy)

Computes and updates the canopy-level conductance (units of m/s) according to the Medlyn model.

The moisture stress factor is applied to `An` already.
"""
function update_canopy_conductance!(p, Y, model::MedlynConductanceModel, canopy)
    c_co2_air = p.drivers.c_co2
    P_air = p.drivers.P
    T_air = p.drivers.T
    q_air = p.drivers.q
    earth_param_set = canopy.parameters.earth_param_set
    thermo_params = earth_param_set.thermo_params
    (; g1, g0, Drel) = canopy.conductance.parameters
    area_index = p.canopy.hydraulics.area_index
    LAI = area_index.leaf
    An = p.canopy.photosynthesis.An
    R = LP.gas_constant(earth_param_set)

    medlyn_factor = @. lazy(medlyn_term(g1, T_air, P_air, q_air, thermo_params))
    @. p.canopy.conductance.r_stomata_canopy =
        1 / upscale_leaf_conductance(
            medlyn_conductance(g0, Drel, medlyn_factor, An, c_co2_air), #conductance, leaf level
            LAI,
            T_air,
            R,
            P_air,
        )
end

# For interfacing with ClimaParams

"""
    function MedlynConductanceParameters(FT::AbstractFloat;
        g1 = 790,
        kwargs...
    )
    function MedlynConductanceParameters(toml_dict;
        g1 = 790,
        kwargs...
    )

Floating-point and toml dict based constructor supplying default values
for the MedlynConductanceParameters struct.
Additional parameter values can be directly set via kwargs.
"""
MedlynConductanceParameters(::Type{FT}; kwargs...) where {FT <: AbstractFloat} =
    MedlynConductanceParameters(CP.create_toml_dict(FT); kwargs...)

function MedlynConductanceParameters(
    toml_dict::CP.AbstractTOMLDict;
    g1 = 790,
    kwargs...,
)
    name_map = (;
        :relative_diffusivity_of_water_vapor => :Drel,
        :min_stomatal_conductance => :g0,
    )

    parameters = CP.get_parameter_values(toml_dict, name_map, "Land")
    FT = CP.float_type(toml_dict)
    g1 = FT.(g1)
    G1 = typeof(g1)
    return MedlynConductanceParameters{FT, G1}(; g1, parameters..., kwargs...)
end

#################### P model conductance ####################
"""
    PModelConductanceParameters{FT <: AbstractFloat}

The required parameters for the Medlyn stomatal conductance model.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct PModelConductanceParameters{
    FT <: AbstractFloat,
}
    "Relative diffusivity of water vapor (unitless)"
    Drel::FT
end

Base.eltype(::PModelConductanceParameters{FT}) where {FT} = FT

struct PModelConductanceParameters{FT, PMCP <: PModelConductanceParameters{FT}} <:
       AbstractStomatalConductanceModel{FT}
    parameters::PMCP
end

function PModelConductance{FT}(
    parameters::PModelConductanceParameters{FT},
) where {FT <: AbstractFloat}
    return PModelConductance{eltype(parameters), typeof(parameters)}(
        parameters,
    )
end

ClimaLand.auxiliary_vars(model::PModelConductance) = (:r_stomata_canopy,)
ClimaLand.auxiliary_types(model::PModelConductance{FT}) where {FT} = (FT,)
ClimaLand.auxiliary_domain_names(::PModelConductance) = (:surface,)

"""
    update_canopy_conductance!(p, Y, model::PModelConductance, canopy)

Computes and updates the canopy-level conductance (units of m/s) according to the P model. 

The moisture stress factor is applied to `An` already.
"""
function update_canopy_conductance!(p, Y, model::PModelConductance, canopy)
    c_co2_air = p.drivers.c_co2
    P_air = p.drivers.P
    T_air = p.drivers.T
    earth_param_set = canopy.parameters.earth_param_set
    (; Drel) = canopy.conductance.parameters
    area_index = p.canopy.hydraulics.area_index
    LAI = area_index.leaf
    ci = p.canopy.photosynthesis.IntVars.ci
    An = p.canopy.photosynthesis.An
    R = LP.gas_constant(earth_param_set)

    @. p.canopy.conductance.r_stomata_canopy =
        1 / upscale_leaf_conductance(
            pmodel_gs(ci / c_co2_air, 
                c_co2_air, 
                An) * Drel, # leaf level conductance in mol H2O Pa^-1
            LAI,
            T_air,
            R,
            P_air,
        )

end

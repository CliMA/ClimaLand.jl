export MedlynConductanceParameters,
    MedlynConductanceModel,
    PModelConductanceParameters,
    PModelConductance,
    CorrectedConductance

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

The moisture stress factor is applied to `An_leaf` already.
"""
function update_canopy_conductance!(p, Y, model::MedlynConductanceModel, canopy)
    c_co2_air = p.drivers.c_co2
    P_air = p.drivers.P
    T_air = p.drivers.T
    q_air = p.drivers.q
    earth_param_set = canopy.earth_param_set
    thermo_params = earth_param_set.thermo_params
    (; g1, g0, Drel) = canopy.conductance.parameters
    area_index = p.canopy.biomass.area_index
    LAI = area_index.leaf
    An_leaf = get_An_leaf(p, canopy.photosynthesis)
    R = LP.gas_constant(earth_param_set)
    FT = typeof(R)
    medlyn_factor = @. lazy(medlyn_term(g1, T_air, P_air, q_air, thermo_params))
    @. p.canopy.conductance.r_stomata_canopy =
        1 / (
            conductance_molar_flux_to_m_per_s(
                medlyn_conductance(g0, Drel, medlyn_factor, An_leaf, c_co2_air), #conductance, leaf level
                T_air,
                R,
                P_air,
            ) * max(LAI, sqrt(eps(FT)))
        ) # multiply by LAI treating all leaves as if they are in parallel
end

# For interfacing with ClimaParams

"""
    function MedlynConductanceParameters(
        toml_dict::CP.ParamDict;
        g1,
        g0 = toml_dict["min_stomatal_conductance"],
    )

TOML dict based constructor supplying default values for the
`MedlynConductanceParameters` struct.
"""
function MedlynConductanceParameters(
    toml_dict::CP.ParamDict;
    g1,
    g0 = toml_dict["min_stomatal_conductance"],
)
    name_map = (; :relative_diffusivity_of_water_vapor => :Drel,)

    parameters = CP.get_parameter_values(toml_dict, name_map, "Land")
    FT = CP.float_type(toml_dict)
    g1 = FT.(g1)
    G1 = typeof(g1)
    return MedlynConductanceParameters{FT, G1}(; g0, g1, parameters...)
end


#################### P model conductance ####################
"""
    PModelConductanceParameters{FT <: AbstractFloat}

The required parameters for the P-Model stomatal conductance model.
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct PModelConductanceParameters{FT <: AbstractFloat}
    "Relative diffusivity of water vapor (unitless)"
    Drel::FT
end

Base.eltype(::PModelConductanceParameters{FT}) where {FT} = FT

struct PModelConductance{FT, PMCP <: PModelConductanceParameters{FT}} <:
       AbstractStomatalConductanceModel{FT}
    parameters::PMCP
end

function PModelConductance{FT}(
    parameters::PModelConductanceParameters{FT},
) where {FT <: AbstractFloat}
    return PModelConductance{eltype(parameters), typeof(parameters)}(parameters)
end

ClimaLand.auxiliary_vars(model::PModelConductance) = (:r_stomata_canopy,)
ClimaLand.auxiliary_types(model::PModelConductance{FT}) where {FT} = (FT,)
ClimaLand.auxiliary_domain_names(::PModelConductance) = (:surface,)

"""
    update_canopy_conductance!(p, Y, model::PModelConductance, canopy)

Computes and updates the canopy-level conductance (units of m/s) according to the P model. 
The P-model predicts the ratio of plant internal to external CO2 concentration χ, and therefore
the stomatal conductance can be inferred from their difference and the net assimilation rate `An`. 

Note that the moisture stress factor `βm` is applied instantaneously to `An` and `gs_co2` in the
P-model photosynthesis update, so it is not applied again here.
"""
function update_canopy_conductance!(p, Y, model::PModelConductance, canopy)
    P_air = p.drivers.P
    T_air = p.drivers.T
    earth_param_set = canopy.earth_param_set
    (; Drel) = canopy.conductance.parameters
    R = LP.gas_constant(earth_param_set)
    FT = eltype(model.parameters)
    @. p.canopy.conductance.r_stomata_canopy =
        1 / (
            conductance_molar_flux_to_m_per_s(
                Drel * p.canopy.photosynthesis.instantaneous.gs_co2, # canopy level conductance in mol H2O/m^2/s
                T_air,
                R,
                P_air,
            ) + eps(FT)
        ) # avoids division by zero, since conductance is zero when An is zero 
end

#################### Corrected conductance ####################
"""
    CorrectedConductance{FT, M, C} <: AbstractStomatalConductanceModel{FT}

A stomatal conductance model `model` whose canopy conductance is multiplied by
the bounded factor `correction`, a [`LogLinearFactor`](@ref) of the canopy
correction features [`canopy_correction_features`](@ref) (leaf area index,
cosine of the solar zenith angle, vapor pressure deficit, snow cover fraction,
top-layer soil water content, moisture stress factor, and log canopy height).

The factor carries empirical corrections of the canopy water-use efficiency,
for example a regression of the evaporative-fraction residuals of the model
against flux-tower observations, into the model without changing its
structure: the corrected conductance enters the same canopy energy and water
balance, so both remain closed.
$(DocStringExtensions.FIELDS)
"""
struct CorrectedConductance{
    FT,
    M <: AbstractStomatalConductanceModel{FT},
    C <: LogLinearFactor{FT},
} <: AbstractStomatalConductanceModel{FT}
    "The stomatal conductance model being corrected"
    model::M
    "The multiplicative correction of the canopy conductance"
    correction::C
end

CorrectedConductance(
    model::AbstractStomatalConductanceModel{FT},
    correction,
) where {FT} = CorrectedConductance{FT, typeof(model), typeof(correction)}(
    model,
    correction,
)

# The wrapped model's `parameters` are reached through the wrapper, as the
# update functions read `canopy.conductance.parameters`.
Base.getproperty(m::CorrectedConductance, s::Symbol) =
    s === :parameters ? getfield(m, :model).parameters : getfield(m, s)

ClimaLand.auxiliary_vars(m::CorrectedConductance) =
    ClimaLand.auxiliary_vars(m.model)
ClimaLand.auxiliary_types(m::CorrectedConductance) =
    ClimaLand.auxiliary_types(m.model)
ClimaLand.auxiliary_domain_names(m::CorrectedConductance) =
    ClimaLand.auxiliary_domain_names(m.model)

"""
    update_canopy_conductance!(p, Y, model::CorrectedConductance, canopy)

Updates the canopy conductance with the wrapped model and divides the canopy
stomatal resistance by the correction factor.
"""
function update_canopy_conductance!(p, Y, model::CorrectedConductance, canopy)
    update_canopy_conductance!(p, Y, model.model, canopy)
    f = model.correction
    x = canopy_correction_features(p, canopy)
    @. p.canopy.conductance.r_stomata_canopy /= f(
        x.LAI,
        x.cosθs,
        x.VPD,
        x.snow_cover_fraction,
        x.θ_top,
        x.βm,
        x.log_height,
    )
end

"""
    canopy_correction_features(p, canopy)

Return a NamedTuple of (lazy) fields of the state variables that the
[`LogLinearFactor`](@ref) corrections of the canopy are functions of: the leaf
area index `LAI`, the cosine of the solar zenith angle `cosθs`, the vapor
pressure deficit `VPD` (kPa), the `snow_cover_fraction` (zero without a snow
model), the volumetric water content of the top soil layer `θ_top`, the
moisture stress factor `βm`, and the logarithm of the canopy height
`log_height` (m, floored at 0.05 m).
"""
function canopy_correction_features(p, canopy)
    thermo_params = LP.thermodynamic_parameters(canopy.earth_param_set)
    LAI = p.canopy.biomass.area_index.leaf
    cosθs = p.drivers.cosθs
    VPD = @. lazy(
        Thermodynamics.vapor_pressure_deficit(
            thermo_params,
            p.drivers.T,
            p.drivers.P,
            p.drivers.q,
        ) / 1000,
    )
    components = Val(canopy.boundary_conditions.prognostic_land_components)
    snow_cover_fraction = canopy_snow_cover_fraction(p, components, LAI)
    θ_top = canopy_top_soil_water(p, canopy.boundary_conditions.ground, LAI)
    βm = p.canopy.soil_moisture_stress.βm
    log_height = log_canopy_height(canopy.biomass.height, eltype(LAI))
    return (; LAI, cosθs, VPD, snow_cover_fraction, θ_top, βm, log_height)
end

"""
    log_canopy_height(height, FT)

The logarithm of the canopy `height` (a number or a field), floored at 0.05 m.
"""
log_canopy_height(height::Number, FT) = log(max(FT(height), FT(0.05)))
log_canopy_height(height, FT) = @. lazy(log(max(height, FT(0.05))))

"""
    canopy_snow_cover_fraction(p, components::Val, LAI)

The snow cover fraction seen by the canopy: that of the snow model when the
land components include `:snow`, zero otherwise.
"""
canopy_snow_cover_fraction(p, ::Val{components}, LAI) where {components} =
    :snow in components ? p.snow.snow_cover_fraction : @. lazy(zero(LAI))

"""
    canopy_top_soil_water(p, ground, LAI)

The volumetric water content of the top soil layer seen by the canopy: that
of the soil model with `PrognosticGroundConditions`, the prescribed driver
`p.drivers.θ` with `PrescribedGroundConditions`.
"""
canopy_top_soil_water(p, ::ClimaLand.PrognosticGroundConditions, LAI) =
    ClimaLand.Domains.top_center_to_surface(p.soil.θ_l)
canopy_top_soil_water(p, ::PrescribedGroundConditions, LAI) = p.drivers.θ

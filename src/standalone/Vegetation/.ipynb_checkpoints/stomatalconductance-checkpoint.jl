export MedlynConductanceParameters,
    MedlynConductanceModel, PModelConductanceParameters, PModelConductance, OWUSConductanceParameters, OWUSConductanceModel

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
    earth_param_set = canopy.parameters.earth_param_set
    thermo_params = earth_param_set.thermo_params
    (; g1, g0, Drel) = canopy.conductance.parameters
    area_index = p.canopy.hydraulics.area_index
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
    function MedlynConductanceParameters(::Type{FT};
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
    toml_dict::CP.ParamDict;
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

Note that the moisture stress factor `βm` is applied to `An` already, so it is not applied again here. 
"""
function update_canopy_conductance!(p, Y, model::PModelConductance, canopy)
    c_co2_air = p.drivers.c_co2
    P_air = p.drivers.P
    T_air = p.drivers.T
    earth_param_set = canopy.parameters.earth_param_set
    (; Drel) = canopy.conductance.parameters
    area_index = p.canopy.hydraulics.area_index
    LAI = area_index.leaf
    ci = p.canopy.photosynthesis.ci             # internal CO2 partial pressure, Pa 
    An_canopy = p.canopy.photosynthesis.An          # net assimilation rate, mol m^-2 s^-1, canopy level
    R = LP.gas_constant(earth_param_set)
    FT = eltype(model.parameters)

    χ = @. lazy(ci / (c_co2_air * P_air))       # ratio of intercellular to ambient CO2 concentration, unitless
    @. p.canopy.conductance.r_stomata_canopy =
        1 / (
            conductance_molar_flux_to_m_per_s(
                gs_h2o_pmodel(χ, c_co2_air, An_canopy, Drel), # canopy level conductance in mol H2O/m^2/s
                T_air,
                R,
                P_air,
            ) + eps(FT)
        ) # avoids division by zero, since conductance is zero when An is zero 
end

#################### OWUS conductance ####################
Base.@kwdef struct OWUSConductanceParameters{
    FT <: AbstractFloat,
}
    "Well-watered transpiration fraction (unitless), plateau of β(s)"
    fww::FT
    "Soil saturation threshold where down-regulation begins (unitless)"
    s_star::FT
    "Soil saturation at shutdown (unitless)"
    s_w::FT
    "Optional cap on stomatal conductance (mol m^-2 s^-1); Inf for none"
    gsw_max::FT = FT(Inf)
end
Base.eltype(::OWUSConductanceParameters{FT}) where {FT} = FT

struct OWUSConductanceModel{FT, OCP <: OWUSConductanceParameters{FT}} <:
       AbstractStomatalConductanceModel{FT}
    parameters::OCP
end

function OWUSConductanceModel{FT}(
    parameters::OWUSConductanceParameters{FT},
) where {FT <: AbstractFloat}
    return OWUSConductanceModel{eltype(parameters), typeof(parameters)}(parameters)
end

"""
    OWUSConductanceModel{FT}(
        ; canopy_params, soil_params, root_params, met_params=nothing, gsw_max = FT(Inf)
    )

Diagnose (fww, s*, s_w) from ClimaLand-style parameter structs and build a model.
"""
function OWUSConductanceModel{FT}(;
    canopy_params, soil_params, root_params, met_params=nothing, gsw_max = FT(Inf)
) where {FT <: AbstractFloat}
    owus = OWUSStomata.build_owus_from_ClimaLand(;
        canopy_params = canopy_params,
        soil_params   = soil_params,
        root_params   = root_params,
        met_params    = met_params,
        overrides     = NamedTuple()
    )
    pars = OWUSConductanceParameters{FT}(;
        fww     = FT(owus.fww),
        s_star  = FT(owus.s_star),
        s_w     = FT(owus.s_w),
        gsw_max = FT(owus.gsw_max),
    )
    return OWUSConductanceModel{FT}(pars)
end


ClimaLand.auxiliary_vars(::OWUSConductanceModel) = (:r_stomata_canopy,)
ClimaLand.auxiliary_types(::OWUSConductanceModel{FT}) where {FT} = (FT,)
ClimaLand.auxiliary_domain_names(::OWUSConductanceModel) = (:surface,)


# --- soft getters: try multiple places, fall back gracefully ---

@inline _has(x, f::Symbol) = Base.hasproperty(x, f)
@inline _get(x, f::Symbol, dflt=nothing) = _has(x, f) ? getproperty(x,f) : dflt

# Soil saturation s \in [0,1].
# Tries: p.soil.θ with soil params (ν, θ_r); else piecewise-stress cache; else 1.0
function _get_saturation(p, canopy)::Float64
    # try full soil state
    soil = _get(p, :soil)
    if soil !== nothing
        θ = _get(soil, :θ)
        params = _get(soil, :parameters)
        if θ !== nothing && params !== nothing
            ν   = _get(params, :ν)
            θ_r = _get(params, :θ_r)
            if ν !== nothing && θ_r !== nothing
                return clamp((mean(ClimaCore.Fields.field_values(θ)) - θ_r) / max(ν - θ_r, eps()), 0.0, 1.0)
            end
        end
    end
    # try canopy stress model cache (θ_low/θ_high)
    sms = _get(p.canopy, :soil_moisture_stress, nothing)
    if sms !== nothing
        θ    = _get(sms, :θ)
        θ_hi = _get(sms, :θ_high)
        θ_lo = _get(sms, :θ_low)
        if θ !== nothing && θ_hi !== nothing && θ_lo !== nothing
            θ̄ = mean(ClimaCore.Fields.field_values(θ))
            return clamp((θ̄ - θ_lo) / max(θ_hi - θ_lo, eps()), 0.0, 1.0)
        end
    end
    return 1.0 # well-watered fallback
end

# Potential evaporation E0 (m s^-1): try energy cache, drivers, or a mild default
function _get_E0_mps(p, canopy)::Float64
    # energy component might expose PET/E0
    energy = canopy.energy
    if _has(p.canopy, :energy)
        E0f = _get(p.canopy.energy, :E0, nothing) # already m s^-1?
        E0d = _get(p.canopy.energy, :E0_mps, nothing)
        if E0f !== nothing
            return Float64(E0f)
        elseif E0d !== nothing
            return Float64(E0d)
        end
    end
    # drivers sometimes carry it
    E0d = _get(p.drivers, :E0, nothing)
    if E0d !== nothing
        return Float64(E0d)
    end
    # gentle default ~2.5 mm/day
    return 2.5e-3 / 86400.0
end

# VPD from (T, P, q); fall back to cached VPD if present
@inline function _svp_pa(Tk::Float64)
    # Tetens (Pa)
    Tc = Tk - 273.15
    return 610.94 * exp((17.625 * Tc) / (Tc + 243.04))
end
function _get_VPD_pa(p)::Float64
    VPD = _get(p.drivers, :VPD, nothing)
    if VPD !== nothing
        return Float64(VPD)
    end
    T = Float64(_get(p.drivers, :T, 298.15))
    P = Float64(_get(p.drivers, :P, 101325.0))
    q = Float64(_get(p.drivers, :q, 0.010))  # kg/kg
    ε = 0.622
    e = (q * P) / (ε + (1 - ε) * q)         # vapor pressure, Pa
    es = _svp_pa(T)
    return max(es - e, 0.0)
end


"""
    update_canopy_conductance!(p, Y, model::OWUSConductanceModel, canopy)

OWUS (Bassiouni–Manzoni–Vico 2023) conductance at canopy scale.
- Diagnose β(s) from params (fww, s*, s_w) and soil saturation `s` from caches.
- Target T = β(s) * E0; invert Fick's law with current VPD,P to get g_sw (mol m^-2 s^-1).
- Convert to m s^-1 and aggregate leaves in parallel using LAI, then write r_stomata_canopy.
"""
function update_canopy_conductance!(p, Y, model::OWUSConductanceModel, canopy)
    FT = eltype(model.parameters)
    earth_param_set = canopy.parameters.earth_param_set
    R = LP.gas_constant(earth_param_set)

    # Drivers and canopy scalars
    P_air = Float64(p.drivers.P)     # Pa
    T_air = Float64(p.drivers.T)     # K
    area_index = p.canopy.hydraulics.area_index
    LAI = Float64(area_index.leaf)

    # Hydromet diagnostics
    s   = _get_saturation(p, canopy)         # 0..1
    E0  = _get_E0_mps(p, canopy)             # m s^-1
    VPD = _get_VPD_pa(p)                     # Pa

    # Use the already-included OWUSStomata module for the leaf-level g_sw
    pars = model.parameters
    owus = OWUSStomata.OWUSStomatalModel(; fww=pars.fww, s_star=pars.s_star, s_w=pars.s_w, gsw_max=pars.gsw_max)
    gsw_ref = Ref{Float64}(0.0)
    OWUSStomata.stomatal_conductance!(gsw_ref, owus;
        s = s, E0 = E0, VPD = VPD, P_air = P_air, T_air = T_air, leaf_factor = 1.0)

    # convert leaf-level g_sw [mol m^-2 s^-1] -> [m s^-1] using ideal gas, then scale to canopy
    g_leaf_mps = gsw_ref[] * (Float64(R) * T_air / P_air)
    FTg = FT( g_leaf_mps * max(LAI, sqrt(eps(Float64))) ) # leaves in parallel
    # write canopy resistance
    @. p.canopy.conductance.r_stomata_canopy = 1 / (FTg + eps(FT))
    return nothing
end

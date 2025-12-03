export MedlynConductanceParameters,
    MedlynConductanceModel, PModelConductanceParameters, PModelConductance, uSPACConductanceParameters, uSPACConductanceModel, uSPACCWDStaticParameters, uSPACConductanceCWDStatic, uSPACPiParameters, uSPACConductancePi

abstract type AbstractStomatalConductanceModel{FT} <:
              AbstractCanopyComponent{FT} end

"""x
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

Note that the moisture stress factor `βm` is applied to `An` already, so it is not applied again here. 
"""
function update_canopy_conductance!(p, Y, model::PModelConductance, canopy)
    c_co2_air = p.drivers.c_co2
    P_air = p.drivers.P
    T_air = p.drivers.T
    earth_param_set = canopy.parameters.earth_param_set
    (; Drel) = canopy.conductance.parameters
    area_index = p.canopy.biomass.area_index
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

#################### uSPAC conductance (existing static back-compat) ####################
Base.@kwdef struct uSPACConductanceParameters{FT <: AbstractFloat}
    "Well-watered transpiration fraction (unitless), plateau of β(s)"
    fww::FT
    "Soil saturation threshold where down-regulation begins (unitless)"
    s_star::FT
    "Soil saturation at shutdown (unitless)"
    s_w::FT
    "Optional cap on stomatal conductance (mol m^-2 s^-1); Inf for none"
    gsw_max::FT = FT(Inf)
end
Base.eltype(::uSPACConductanceParameters{FT}) where {FT} = FT

struct uSPACConductanceModel{FT, OCP <: uSPACConductanceParameters{FT}} <:
       AbstractStomatalConductanceModel{FT}
    parameters::OCP
end

function uSPACConductanceModel{FT}(
    parameters::uSPACConductanceParameters{FT},
) where {FT <: AbstractFloat}
    return uSPACConductanceModel{eltype(parameters), typeof(parameters)}(parameters)
end

function uSPACConductanceModel{FT}(;
    canopy_params, soil_params, root_params, met_params=nothing, gsw_max = FT(Inf)
) where {FT <: AbstractFloat}
    uspac = uSPACStomata.build_uspac_from_ClimaLand(;
        canopy_params = canopy_params,
        soil_params   = soil_params,
        root_params   = root_params,
        met_params    = met_params,
        overrides     = NamedTuple()
    )
    pars = uSPACConductanceParameters{FT}(;
        fww     = FT(uspac.fww),
        s_star  = FT(uspac.s_star),
        s_w     = FT(uspac.s_w),
        gsw_max = FT(uspac.gsw_max),
    )
    return uSPACConductanceModel{FT}(pars)
end


ClimaLand.auxiliary_vars(::uSPACConductanceModel) = (:r_stomata_canopy, :gsw_leaf, :gsw_canopy)
ClimaLand.auxiliary_types(::uSPACConductanceModel{FT}) where {FT} = (FT, FT, FT)
ClimaLand.auxiliary_domain_names(::uSPACConductanceModel) = (:surface, :surface, :surface)

# --- soft getters ---
@inline _has(x, f::Symbol) = Base.hasproperty(x, f)
@inline _get(x, f::Symbol, dflt=nothing) = _has(x, f) ? getproperty(x,f) : dflt

# Soil saturation s \in [0,1].
# Field-valued soil saturation s ∈ [0,1]
@inline function _saturation_field(p, ::Any, ::Type{FT}) where {FT}
    # pick any canopy Field to clone shape/grid from
    like = p.drivers.T

    # 1) Prefer canopy soil_moisture_stress if it exposes θ info
    if Base.hasproperty(p.canopy, :soil_moisture_stress)
        sms = p.canopy.soil_moisture_stress

        if Base.hasproperty(sms, :θ) &&
           Base.hasproperty(sms, :θ_high) &&
           Base.hasproperty(sms, :θ_low)
            θ    = sms.θ
            θ_hi = sms.θ_high
            θ_lo = sms.θ_low
            return @. clamp((θ - θ_lo) / max(θ_hi - θ_lo, eps(FT)), FT(0), FT(1))
        elseif Base.hasproperty(sms, :s)
            s = sms.s
            return @. clamp(FT(s), FT(0), FT(1))
        elseif Base.hasproperty(sms, :βm)
            βm = sms.βm
            # treat βm as a saturation-like 0..1 proxy
            return @. clamp(FT(βm), FT(0), FT(1))
        end
    end

    # 2) Otherwise, try soil state with θ_r, ν
    if Base.hasproperty(p, :soil)
        soil = p.soil
        if Base.hasproperty(soil, :θ) && Base.hasproperty(soil, :parameters)
            θ = soil.θ
            pars = soil.parameters
            if Base.hasproperty(pars, :ν) && Base.hasproperty(pars, :θ_r)
                ν   = pars.ν
                θ_r = pars.θ_r
                return @. clamp((θ - θ_r) / max(ν - θ_r, eps(FT)), FT(0), FT(1))
            end
        end
    end
end

@inline function _E0_field(p, ::Any, ::Type{FT}) where {FT}
    if Base.hasproperty(p.canopy, :energy) && Base.hasproperty(p.canopy.energy, :E0)
        p.canopy.energy.E0
    elseif Base.hasproperty(p.drivers, :E0)
        p.drivers.E0
    else
        # fill a Field with a small constant default
        LAI = p.canopy.biomass.area_index.leaf
        @. zero(LAI) + FT(2.5e-3/86400)
    end
end

@inline function _VPD_field(p, ::Type{FT}) where {FT}
    T = p.drivers.T
    P = p.drivers.P
    q = p.drivers.q
    @inline _svp_pa(Tk) = FT(610.94) * exp((FT(17.625) * (Tk - FT(273.15))) / (Tk - FT(273.15) + FT(243.04)))
    e  = @. (q * P) / (FT(0.622) + (FT(1) - FT(0.622)) * q)
    es = @. _svp_pa(T)
    @. max(es - e, FT(0))
end

"""
    update_canopy_conductance!(p, Y, model::uSPACConductanceModel, canopy)

Static uSPAC (back-compat, no time variability in parameters).
"""
function update_canopy_conductance!(p, Y, model::uSPACConductanceModel, canopy)
    FT   = eltype(model.parameters)
    pars = model.parameters
    earth = canopy.parameters.earth_param_set
    Rgas  = LP.gas_constant(earth)  # FT scalar

    # Inputs as Fields
    P_air = p.drivers.P
    T_air = p.drivers.T
    LAI   = p.canopy.biomass.area_index.leaf
    s     = _saturation_field(p, canopy, FT)
    E0    = _E0_field(p, canopy, FT)
    VPD   = _VPD_field(p, FT)

    # uSPAC β(s)
    r  = @. clamp((s - pars.s_w) / max(pars.s_star - pars.s_w, eps(FT)), FT(0), FT(1))
    β  = @. pars.fww * r

    # gsw (leaf molar), all Fields
    ρw = FT(1000.0); Mw = FT(0.01801528)
    E_mps = @. β * E0
    E_mol = @. (ρw * E_mps) / Mw
    gsw_leaf = @. ifelse(VPD > eps(FT), min(E_mol * (P_air / VPD), pars.gsw_max), FT(0))

    # expose leaf/canopy molar conductance (optional aux)
    if Base.hasproperty(p.canopy.conductance, :gsw_leaf)
        @. p.canopy.conductance.gsw_leaf = gsw_leaf
    end
    if Base.hasproperty(p.canopy.conductance, :gsw_canopy)
        @. p.canopy.conductance.gsw_canopy = gsw_leaf * max(LAI, sqrt(eps(FT)))
    end

    # convert to m s^-1 and to canopy resistance Field
    g_leaf_mps = @. gsw_leaf * (Rgas * T_air / P_air)             # m s^-1 (leaf)
    g_canopy   = @. g_leaf_mps * max(LAI, sqrt(eps(FT)))           # m s^-1 (ground)
    @. p.canopy.conductance.r_stomata_canopy = 1 / (g_canopy + eps(FT))
    return nothing
end


# ===================== CWD-controlled uSPAC =====================

# Maps w = [1, log1p(CWD_mm)] -> (α, a, b) via Γ (rows α,a,b), then
#   fww = σ(α), s_w = σ(a), s* = s_w + (1 - s_w) σ(b)
Base.@kwdef struct uSPACCWDStaticParameters{FT <: AbstractFloat}
    Γ::NTuple{6,FT}       # row-major [α0, αCWD, a0, aCWD, b0, bCWD]
    cwd_mm::FT            # site/pixel covariate (mm)
    gsw_max::FT = FT(Inf)
end
Base.eltype(::uSPACCWDStaticParameters{FT}) where {FT} = FT

struct uSPACConductanceCWDStatic{FT,
                                P<:uSPACCWDStaticParameters{FT}} <:
       AbstractStomatalConductanceModel{FT}
    parameters::P
end
uSPACConductanceCWDStatic{FT}(p::uSPACCWDStaticParameters{FT}) where {FT<:AbstractFloat} =
    uSPACConductanceCWDStatic{FT,typeof(p)}(p)


ClimaLand.auxiliary_vars(::uSPACConductanceCWDStatic) = (:r_stomata_canopy, :gsw_leaf, :gsw_canopy)
ClimaLand.auxiliary_types(::uSPACConductanceCWDStatic{FT}) where {FT} = (FT, FT, FT)
ClimaLand.auxiliary_domain_names(::uSPACConductanceCWDStatic) = (:surface, :surface, :surface)

@inline _σ(x) = inv(one(x) + exp(-x))
@inline _affine2(γ0::T, γ1::T, CWD::T) where {T} = γ0 + γ1 * log1p(CWD)

function update_canopy_conductance!(p, Y, model::uSPACConductanceCWDStatic, canopy)
    FT   = eltype(model.parameters)
    pars = model.parameters
    earth = canopy.parameters.earth_param_set
    Rgas  = LP.gas_constant(earth)

    # Drivers / state
    P_air = p.drivers.P
    T_air = p.drivers.T
    LAI   = p.canopy.biomass.area_index.leaf
    s     = _saturation_field(p, canopy, FT)
    E0    = _E0_field(p, canopy, FT)
    VPD   = _VPD_field(p, FT)

    # CWD → (α, a, b) → (fww, s*, sw)  (scalars of type FT)
    γα0, γαC, γa0, γaC, γb0, γbC = pars.Γ
    CWD = FT(pars.cwd_mm)
    α = γα0 + γαC * log1p(CWD)
    a = γa0  + γaC  * log1p(CWD)
    b = γb0  + γbC  * log1p(CWD)

    sw    = inv(FT(1) + exp(-a))
    sb    = inv(FT(1) + exp(-b))
    sstar = sw + (FT(1) - sw) * sb
    fww   = inv(FT(1) + exp(-α))

    # β(s) piecewise (broadcasted)
    r  = @. clamp((s - sw) / max(sstar - sw, eps(FT)), FT(0), FT(1))
    β  = @. fww * r

    # Convert to gsw (leaf, molar)
    ρw = FT(1000.0); Mw = FT(0.01801528)
    E_mps = @. β * E0
    E_mol = @. (ρw * E_mps) / Mw
    gsw_max = pars.gsw_max
    gsw_leaf = @. ifelse(VPD > eps(FT),
                         min(E_mol * (P_air / VPD), gsw_max),
                         FT(0))

    if Base.hasproperty(p.canopy.conductance, :gsw_leaf)
        @. p.canopy.conductance.gsw_leaf = gsw_leaf
    end
    if Base.hasproperty(p.canopy.conductance, :gsw_canopy)
        @. p.canopy.conductance.gsw_canopy = gsw_leaf * max(LAI, sqrt(eps(FT)))
    end

    g_leaf_mps = @. gsw_leaf * (Rgas * T_air / P_air)
    g_canopy   = @. g_leaf_mps * max(LAI, sqrt(eps(FT)))
    @. p.canopy.conductance.r_stomata_canopy = 1 / (g_canopy + eps(FT))
    return nothing
end

# ===================== Π-group controlled uSPAC =====================

Base.@kwdef struct uSPACPiParameters{FT <: AbstractFloat}
    # Coefficients to map aridity → ΠR, ΠF  (Π = α + β * log1p(aridity))
    ΓR::NTuple{2,FT}   # [α_R, β_R]
    ΓF::NTuple{2,FT}   # [α_F, β_F]

    # Coefficients to map soil texture → ΠT, ΠS  (Π = α + β*Sand)
    ΓT::NTuple{2,FT}   # sand    → ΠT  
    ΓS::NTuple{2,FT}   # sand    → ΠS

    # Covariates can be scalars (FT) or Field/array-like; we keep them as Any to allow either.
    #aridity_idx::Any   # e.g., CWD (mm), PET/MAP ratio, etc.  scalar or Field
    #sand::Any          # 0..1 scalar or Field
    aridity_idx::Any = nothing   # <- default
    sand::Any        = nothing   # <- default

    # Soil water retention exponent; thresholds for s* and s_w as fractions of fww
    b::FT = FT(4.38)
    β_star_frac::FT = FT(0.95)
    β_w_frac::FT    = FT(0.05)

    gsw_max::FT = FT(Inf)
end
Base.eltype(::uSPACPiParameters{FT}) where {FT} = FT

struct uSPACConductancePi{FT,
                          P<:uSPACPiParameters{FT}} <:
       AbstractStomatalConductanceModel{FT}
    parameters::P
end
uSPACConductancePi{FT}(p::uSPACPiParameters{FT}) where {FT<:AbstractFloat} =
    uSPACConductancePi{FT,typeof(p)}(p)

ClimaLand.auxiliary_vars(::uSPACConductancePi) = (:r_stomata_canopy, :gsw_leaf, :gsw_canopy)
ClimaLand.auxiliary_types(::uSPACConductancePi{FT}) where {FT} = (FT, FT, FT)
ClimaLand.auxiliary_domain_names(::uSPACConductancePi) = (:surface, :surface, :surface)

# helpers
@inline function _getfirst(x, names::NTuple{N,Symbol}, default=nothing) where {N}
    x === nothing && return default
    for nm in names
        if Base.hasproperty(x, nm)
            return getproperty(x, nm)
        end
    end
    return default
end

# --- helpers to support scalar OR Field covariates ---
@inline _as_like_field(like, x::Number, ::Type{FT}) where {FT} = @. zero(like) + FT(x)
@inline _as_like_field(like, x,       ::Type)               = x   # assume array/Field-like already shaped

#@inline _σ(x) = inv(one(x) + exp(-x)) DEFINED ABOVE ELSEWHERE

# Elementwise uSPAC Π → (fww, s*, s_w); written as scalar funcs so we can broadcast them
# --- covariate extractors (auto-pull; return Fields matching grid) ---
function extract_aridity!(like, p, canopy, ::Type{FT}) where {FT}
    # drivers first; then canopy.parameters; else soft default
    ari = _getfirst(p.drivers, (:CWD_mm, :aridity_index, :PET_MAP_ratio, :MAP_over_PET), nothing)
    return _as_like_field(like, ari, FT)
end

function extract_sand!(like, p, ::Type{FT}) where {FT}
    # try soil, then soil.parameters, then parameters.texture; else default ~loamy
    sand = _getfirst(p.soil, (:ν_ss_quartz, :sand, :soil_sand, :sand_frac), nothing)
    if (sand === nothing) && Base.hasproperty(p.soil, :parameters)
        sand = _getfirst(p.soil.parameters, (:ν_ss_quartz, :sand, :soil_sand, :sand_frac), sand)
        tex  = _getfirst(p.soil.parameters, (:texture,), nothing)
        if tex !== nothing
            sand = sand === nothing ? _getfirst(tex, (:ν_ss_quartz, :sand, :sand_frac), sand) : sand
        end
    end
    sand === nothing && (sand = FT(0.45))  # safe default fraction (0..1)
    return _as_like_field(like, sand, FT)
end

# --- uSPAC algebra as scalar funcs so we can broadcast ---
@inline function _uspac_fww_from_Pi(ΠR, ΠF)
    ΠR_safe = max(abs(ΠR), sqrt(eps(ΠR)))  # Prevent division by zero
    halfΠF = ΠF / 2
    rad = (halfΠF + 1)^2 - 2 * ΠF * ΠR_safe
    rad = max(rad, eps(rad))  # Ensure positive radicand
    fww = 1 - (1 / (2 * ΠR_safe)) * (1 + halfΠF - sqrt(rad))
    return clamp(fww, zero(fww), one(fww))  # Ensure [0,1]
end

@inline function _uspac_s_of_beta(β, ΠR, ΠF, ΠT, ΠS, b)
    β = clamp(β, zero(β), one(β))
    
    # Bassiouni et al Eq. 5
    denom = 1 - (1 - β) * ΠR
    denom = ifelse(abs(denom) < sqrt(eps(denom)), sign(denom) * sqrt(eps(denom)), denom)
    
    # Safeguard ΠT and ΠS
    ΠT_safe = max(abs(ΠT), sqrt(eps(ΠT)))
    ΠS_safe = max(abs(ΠS), sqrt(eps(ΠS)))
    
    termA = (4 * β * ΠS_safe * ΠS_safe) / ΠT_safe
    termB = (2 * (1 - β) - β * ΠF) / denom
    inner = sqrt(max(1 + termA * termB, eps(termA))) - 1
    base  = (ΠT_safe / (2 * β * ΠS_safe)) * inner
    s = (max(base, eps(base)))^(-one(b)/b)
    return clamp(s, zero(s), one(s))
end

function update_canopy_conductance!(p, Y, model::uSPACConductancePi, canopy)
    FT   = eltype(model.parameters)
    pars = model.parameters
    earth = canopy.parameters.earth_param_set
    Rgas  = LP.gas_constant(earth)

    # State/driver Fields
    P_air = p.drivers.P
    T_air = p.drivers.T
    LAI   = p.canopy.biomass.area_index.leaf
    s     = _saturation_field(p, canopy, FT)
    E0    = _E0_field(p, canopy, FT)
    VPD   = _VPD_field(p, FT)

    like = s  # grid shape for broadcasting

    # --- auto-fetch covariates as Fields ---
    # aridity = extract_aridity!(like, p, canopy, FT)  # Field or broadcasted scalar
    # sand    = extract_sand!(like, p, FT)
    like = s
    aridity = isnothing(pars.aridity_idx) ? extract_aridity!(like, p, canopy, FT) : _as_like_field(like, pars.aridity_idx, FT)
    sand    = isnothing(pars.sand) ? extract_sand!(like, p, FT) : _as_like_field(like, pars.sand, FT)

    # --- Π from covariates (sand-only for ΠT, ΠS) ---
    αR, βR  = pars.ΓR
    αF, βF  = pars.ΓF
    αT, βTs = pars.ΓT
    αS, βSs = pars.ΓS

    ΠR = @. αR + βR * log1p(aridity)
    ΠF = @. αF + βF * log1p(aridity)
    ΠT = @. αT + βTs * sand
    ΠS = @. αS + βSs * sand

    # Π → (fww, s*, s_w) elementwise
    fww   = @. _uspac_fww_from_Pi(ΠR, ΠF)
    sstar = @. _uspac_s_of_beta(pars.β_star_frac * fww, ΠR, ΠF, ΠT, ΠS, pars.b)
    sw    = @. _uspac_s_of_beta(pars.β_w_frac    * fww, ΠR, ΠF, ΠT, ΠS, pars.b)

    # β(s) piecewise
    r  = @. clamp((s - sw) / max(sstar - sw, eps(FT)), FT(0), FT(1))
    β  = @. fww * r

    # Stomatal conductance (leaf, molar) and canopy resistance (ground)
    ρw = FT(1000.0); Mw = FT(0.01801528)
    E_mps   = @. β * E0
    E_mol   = @. (ρw * E_mps) / Mw
    gsw_leaf = @. ifelse(VPD > eps(FT),
                         min(E_mol * (P_air / VPD), pars.gsw_max),
                         FT(0))

    if Base.hasproperty(p.canopy.conductance, :gsw_leaf)
        @. p.canopy.conductance.gsw_leaf = gsw_leaf
    end
    if Base.hasproperty(p.canopy.conductance, :gsw_canopy)
        @. p.canopy.conductance.gsw_canopy = gsw_leaf * max(LAI, sqrt(eps(FT)))
    end

    g_leaf_mps = @. gsw_leaf * (Rgas * T_air / P_air)
    g_canopy   = @. g_leaf_mps * max(LAI, sqrt(eps(FT)))
    @. p.canopy.conductance.r_stomata_canopy = 1 / (g_canopy + eps(FT))
    return nothing
end
export MedlynConductanceParameters,
    MedlynConductanceModel, PModelConductanceParameters, PModelConductance, uSPACPiParameters, uSPACConductancePi

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

#################### uSPAC conductance ####################

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

# ===================== Π-group controlled uSPAC =====================

Base.@kwdef struct uSPACPiParameters{FT <: AbstractFloat}
    # Calibrated parameters (aridity-dependent)
    βkx_base::FT      # Log-scale intercept for absolute kx (m day⁻¹ MPa⁻¹)
    βkx_coord::FT     # Coordination coefficient: kx ~ P50^βkx_coord (safety-efficiency trade-off)
    βψx50_base::FT    # Base for ψx50 exponential transform
    βψx50_slope::FT   # Slope for ψx50 vs. aridity
    βΠR_base::FT      # Base for ΠR logistic transform
    βΠR_slope::FT     # Slope for ΠR logistic transform
    
    # Trait variance parameters (Phase 1.5: Climate-dependent variance, best-guess functional forms)
    # Base values from literature (wet climate reference)
    σ_logkx_base::FT = FT(0.5)    # Base std dev in log(kx) (wet climates, Liu et al. 2019)
    σ_P50_base::FT = FT(1.5)      # Base std dev of P50 in MPa (Choat et al. 2012)
    σ_ΠR_base::FT = FT(0.15)      # Base std dev of ΠR (unitless)
    
    # Climate-dependent variance modifiers (aridity → variance relationship)
    # Hypothesis: Dry climates have higher trait variance (niche partitioning)
    α_σ_logkx::FT = FT(0.3)   # Increase in σ_logkx per unit decrease in aridity (dry → more variable kx)
    α_σ_P50::FT = FT(1.0)     # Increase in σ_P50 per unit decrease in aridity (dry → more variable P50)
    α_σ_ΠR::FT = FT(0.1)      # Increase in σ_ΠR per unit decrease in aridity (dry → more variable strategy)
    
    # Correlations (also climate-dependent)
    ρ_kx_P50_base::FT = FT(0.7)   # Base correlation (Manzoni et al. 2013)
    ρ_kx_ΠR_base::FT = FT(-0.3)   # Base correlation
    ρ_P50_ΠR_base::FT = FT(-0.2)  # Base correlation
    α_ρ_kx_P50::FT = FT(-0.2)     # Weaker coordination in dry climates (more independent strategies)
    
    n_quad::Int = 3                      # Quadrature order (3, 5, or 7 points per dimension)
    use_trait_distribution::Bool = false  # Enable trait heterogeneity
    
    # Covariates
    aridity_idx::Any = nothing
    
    # Fixed parameters
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

    # --- Extract aridity (only covariate needed for calibration) ---
    aridity = isnothing(pars.aridity_idx) ? 
              extract_aridity!(like, p, canopy, FT) : 
              _as_like_field(like, pars.aridity_idx, FT)

    # Normalization: For Global Aridity Index (P/ET0):
    # 0 = hyper-arid, 0.03-0.2 = arid, 0.2-0.5 = semi-arid, 0.5-0.65 = dry sub-humid
    # 0.65-1.0 = humid, >1.0 = very humid, >2.0 = extremely humid
    # We INVERT so that: high norm = dry, low norm = wet
    # Clamp to reasonable range (2.0 = very humid reference)
    # NOTE: Ocean points are already masked to NaN in the aridity field by 
    # load_aridity_field() using the topographic land/sea mask
    aridity_norm = @. 1.0 - clamp(aridity / FT(2.0), FT(0.0), FT(1.0))

    # --- Compute kx, ψx50, and ΠR from CALIBRATED β coefficients ---

    #     WET ←─────────── aridity_norm ─────────→ DRY
    # 0.0                                       1.0
    # (aridity_norm: 0 = wet, 1 = dry)

    # ψx50:  ╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲╲  (exponential, always negative, DECREASES with aridity_norm)
    #     -0.5 MPa (wet) → -10 MPa (dry)

    # kx:    ╲╲╲╲╲╲╲╲╲╲╲╲╲╲  (decreases with aridity_norm via coordination with ψx50)
    #     high → low

    # ΠR:    ───╮
    #           ╰──────  (sigmoid S-curve)
    #     0 (isohydric, wet) → 1 (anisohydric, dry)

    # ψx50: Negative water potential using exponential
    # Higher aridity_norm (dry) → more negative P50 (cavitation-resistant)
    # Clamp the exponent to prevent extreme values in hyperarid regions
    # Use ifelse to skip ocean points (aridity_norm = NaN)
    ψx50_exponent = @. ifelse(isnan(aridity_norm), FT(0), 
                              clamp(pars.βψx50_base + pars.βψx50_slope * aridity_norm, FT(-2.0), FT(4.0)))
    ψx50_mean = @. ifelse(isnan(aridity_norm), FT(NaN), -exp(ψx50_exponent))
    
    # Encode the universal trade-off: safety-efficiency spectrum
    # literature shows kx and P50 are coordinated
    # High kx → vulnerable to cavitation → requires less negative P50 (wet climates)
    # Low kx → safer from cavitation → can have very negative P50 (dry climates)
    # kx's direct climate response is weak once you account for P50 (Liu et al. 2019)
    # 
    # COORDINATION EQUATION: kx = exp(βkx_base - βkx_coord * log(-ψx50))
    # Note the NEGATIVE sign: as ψx50 becomes more negative (larger -ψx50),
    # log(-ψx50) increases, so kx DECREASES (with positive βkx_coord)
    # This creates: less negative P50 → high kx (wet), more negative P50 → low kx (dry)
    # Clamp log argument to prevent issues with extreme ψx50 values
    # Also clamp the full exponent to prevent overflow/underflow in exp()
    kx_exponent = @. pars.βkx_base - pars.βkx_coord * log(clamp(-ψx50_mean, FT(0.1), FT(100.0)))
    kx_mean = @. ifelse(isnan(aridity_norm), FT(NaN),
                        exp(clamp(kx_exponent, FT(-20.0), FT(10.0))))

    # ΠR: Logistic sigmoid (0 to 1)
    # Higher aridity_norm (dry) → higher ΠR (anisohydric strategy)
    # Clamp the sigmoid argument to prevent overflow in exp()
    ΠR_logit = @. pars.βΠR_base + pars.βΠR_slope * aridity_norm
    ΠR_mean = @. ifelse(isnan(aridity_norm), FT(NaN),
                        1 / (1 + exp(-clamp(ΠR_logit, FT(-20.0), FT(20.0)))))

    # ============ TRAIT DISTRIBUTION INTEGRATION (Phase 1.5: Climate-dependent variance) ============
    if pars.use_trait_distribution
        # Compute climate-dependent variance parameters
        # Hypothesis: Dry climates → higher variance (niche partitioning for limited water)
        #             Wet climates → lower variance (competitive exclusion, less niche space)
        # aridity_norm: 0 (dry) → 1 (wet), so (1 - aridity_norm) = aridity stress
        
        # Standard deviations increase with aridity stress (dry conditions)
        # Skip ocean points (aridity_norm = NaN) → set variance to NaN
        aridity_stress = @. ifelse(isnan(aridity_norm), FT(NaN), FT(1) - aridity_norm)
        σ_logkx = @. ifelse(isnan(aridity_stress), FT(NaN), pars.σ_logkx_base + pars.α_σ_logkx * aridity_stress)
        σ_P50 = @. ifelse(isnan(aridity_stress), FT(NaN), pars.σ_P50_base + pars.α_σ_P50 * aridity_stress)
        σ_ΠR = @. ifelse(isnan(aridity_stress), FT(NaN), pars.σ_ΠR_base + pars.α_σ_ΠR * aridity_stress)
        
        # Correlations weaken in dry climates (more independent strategies)
        ρ_kx_P50 = @. ifelse(isnan(aridity_stress), FT(NaN), pars.ρ_kx_P50_base + pars.α_ρ_kx_P50 * aridity_stress)
        ρ_kx_ΠR = pars.ρ_kx_ΠR_base  # Keep fixed for now
        ρ_P50_ΠR = pars.ρ_P50_ΠR_base  # Keep fixed for now
        
        # Clamp correlations to valid range [-1, 1]
        ρ_kx_P50 = @. clamp(ρ_kx_P50, FT(-0.99), FT(0.99))
        
        # For simplicity in quadrature generation, use spatial mean variance
        # (Full spatially-varying variance would require per-gridcell quadrature)
        # Important: Filter out NaN values (ocean points) before computing means
        # Note: We replace NaN with 0 for summing, then count non-NaN points
        σ_logkx_clean = @. ifelse(isnan(σ_logkx), FT(0), σ_logkx)
        σ_P50_clean = @. ifelse(isnan(σ_P50), FT(0), σ_P50)
        σ_ΠR_clean = @. ifelse(isnan(σ_ΠR), FT(0), σ_ΠR)
        ρ_kx_P50_clean = @. ifelse(isnan(ρ_kx_P50), FT(0), ρ_kx_P50)
        aridity_valid_count = @. ifelse(isnan(aridity_norm), FT(0), FT(1))
        n_valid = max(sum(aridity_valid_count), FT(1))
        
        σ_logkx_val = sum(σ_logkx_clean) / n_valid
        σ_P50_val = sum(σ_P50_clean) / n_valid
        σ_ΠR_val = sum(σ_ΠR_clean) / n_valid
        ρ_kx_P50_val = sum(ρ_kx_P50_clean) / n_valid
        ρ_kx_ΠR_val = ρ_kx_ΠR
        ρ_P50_ΠR_val = ρ_P50_ΠR
        
        # Generate trait quadrature (3D: log(kx), P50, ΠR)
        # Mean traits come from climate-dependent calibration
        # Important: Avoid log(NaN) at ocean points - use cleaned fields
        μ_logkx_field = @. ifelse(isnan(aridity_norm), FT(0), log(kx_mean))
        μ_P50_field = @. ifelse(isnan(aridity_norm), FT(0), ψx50_mean)
        μ_ΠR_field = @. ifelse(isnan(aridity_norm), FT(0), ΠR_mean)
        
        # Extract representative mean values for quadrature generation
        # Use spatial mean of valid (non-ocean) points only
        μ_logkx_val = sum(μ_logkx_field) / n_valid
        μ_P50_val = sum(μ_P50_field) / n_valid
        μ_ΠR_val = sum(μ_ΠR_field) / n_valid
        
        # Point-wise integration over trait space
        # We need to broadcast the integration over all spatial points
        gsw_leaf = similar(s)  # Initialize output field
        
        # Generate quadrature points using spatially-averaged variance and means
        # (For full spatial heterogeneity, would need per-gridcell quadrature - future optimization)
        quad = uSPACStomata.generate_trait_quadrature(;
            FT = FT,
            n_quad = pars.n_quad,
            μ_logkx = μ_logkx_val,  # Spatial mean
            μ_P50 = μ_P50_val,
            μ_ΠR = μ_ΠR_val,
            σ_logkx = σ_logkx_val,  # Climate-dependent variance (spatial average)
            σ_P50 = σ_P50_val,
            σ_ΠR = σ_ΠR_val,
            ρ_kx_P50 = ρ_kx_P50_val,  # Climate-dependent correlations
            ρ_kx_ΠR = ρ_kx_ΠR_val,
            ρ_P50_ΠR = ρ_P50_ΠR_val
        )
        
        # Integrate over traits: E[gsw] = Σ wᵢ × gsw(traitᵢ)
        @. gsw_leaf = zero(s)  # Initialize
        for i in 1:length(quad.points)
            logkx_sample, P50_sample, ΠR_sample = quad.points[i]
            kx_sample_val = exp(logkx_sample)
            weight = quad.weights[i]
            
            # Compute Π-groups for this trait sample (broadcast over space)
            kx_sample_field = @. zero(kx_mean) + kx_sample_val
            P50_sample_field = @. zero(ψx50_mean) + P50_sample
            ΠR_sample_field = @. zero(ΠR_mean) + ΠR_sample
            
            ΠF_sample, ΠT_sample, ΠS_sample = uSPACStomata.pi_groups_from_calibrated_traits(
                p, canopy, pars, kx_sample_field, P50_sample_field, ΠR_sample_field, FT
            )
            
            # Π → (fww, s*, s_w) for this trait sample
            fww_sample = @. _uspac_fww_from_Pi(ΠR_sample, ΠF_sample)
            sstar_sample = @. _uspac_s_of_beta(
                pars.β_star_frac * fww_sample, 
                ΠR_sample, ΠF_sample, ΠT_sample, ΠS_sample, pars.b
            )
            sw_sample = @. _uspac_s_of_beta(
                pars.β_w_frac * fww_sample,
                ΠR_sample, ΠF_sample, ΠT_sample, ΠS_sample, pars.b
            )
            
            # β(s) piecewise for this trait sample
            r_sample = @. clamp((s - sw_sample) / max(sstar_sample - sw_sample, eps(FT)), FT(0), FT(1))
            β_sample = @. fww_sample * r_sample
            
            # Stomatal conductance for this trait sample
            ρw = FT(1000.0); Mw = FT(0.01801528)
            E_mps_sample = @. β_sample * E0
            E_mol_sample = @. (ρw * E_mps_sample) / Mw
            gsw_sample = @. ifelse(VPD > eps(FT),
                                  min(E_mol_sample * (P_air / VPD), pars.gsw_max),
                                  FT(0))
            
            # Accumulate weighted contribution
            @. gsw_leaf = gsw_leaf + weight * gsw_sample
        end
    else
        # ============ MEAN-FIELD (Original Code) ============
        # Compute remaining Π-groups from model components
        ΠF, ΠT, ΠS = uSPACStomata.pi_groups_from_calibrated_traits(p, canopy, pars, kx_mean, ψx50_mean, ΠR_mean, FT)

        # --- Π → (fww, s*, s_w) elementwise ---
        fww   = @. _uspac_fww_from_Pi(ΠR_mean, ΠF)
        sstar = @. _uspac_s_of_beta(pars.β_star_frac * fww, ΠR_mean, ΠF, ΠT, ΠS, pars.b)
        sw    = @. _uspac_s_of_beta(pars.β_w_frac    * fww, ΠR_mean, ΠF, ΠT, ΠS, pars.b)

        # --- β(s) piecewise ---
        r  = @. clamp((s - sw) / max(sstar - sw, eps(FT)), FT(0), FT(1))
        β  = @. fww * r

        # --- Stomatal conductance (leaf, molar) and canopy resistance ---
        ρw = FT(1000.0); Mw = FT(0.01801528)
        E_mps   = @. β * E0
        E_mol   = @. (ρw * E_mps) / Mw
        gsw_leaf = @. ifelse(VPD > eps(FT),
                             min(E_mol * (P_air / VPD), pars.gsw_max),
                             FT(0))
    end

    # ============ COMMON OUTPUT (Both mean-field and distribution) ============

    # Optional aux storage
    if Base.hasproperty(p.canopy.conductance, :gsw_leaf)
        @. p.canopy.conductance.gsw_leaf = gsw_leaf
    end
    if Base.hasproperty(p.canopy.conductance, :gsw_canopy)
        @. p.canopy.conductance.gsw_canopy = gsw_leaf * max(LAI, sqrt(eps(FT)))
    end

    # Convert to canopy resistance (ground area basis)
    g_leaf_mps = @. gsw_leaf * (Rgas * T_air / P_air)
    g_canopy   = @. g_leaf_mps * max(LAI, sqrt(eps(FT)))
    @. p.canopy.conductance.r_stomata_canopy = 1 / (g_canopy + eps(FT))
    return nothing
end
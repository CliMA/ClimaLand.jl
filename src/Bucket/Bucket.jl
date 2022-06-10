module Bucket
# (1) Add snow - is there a better way to handle timesteps where only
# some of the surface energy is used to melt all of the snow?
# (2) Make atmos functions a function of space (lat lon? x y? coords?)
# (3) Look at allocations & performance...
using UnPack
using DocStringExtensions
using CLIMAParameters: AbstractEarthParameterSet
using CLIMAParameters: Stefan
using CLIMAParameters.Planet:
    LH_v0, LH_s0, LH_f0, ρ_cloud_liq, T_freeze, T_0, R_d, cp_d, grav, cp_v
using SurfaceFluxes
using SurfaceFluxes.UniversalFunctions
using Thermodynamics
using StaticArrays
using ClimaCore.Fields: coordinate_field
using ClimaLSM
import ClimaLSM.Domains: coordinates, Plane
import ClimaLSM:
    AbstractModel,
    make_update_aux,
    make_rhs,
    prognostic_vars,
    auxiliary_vars,
    name,
    prognostic_types,
    auxiliary_types
export BucketModelParameters,
    PrescribedAtmosphere,
    PrescribedRadiation,
    BucketModel,
    AbstractAtmosphericDrivers,
    AbstractRadiativeDrivers,
    surface_fluxes,
    surface_air_density,
    liquid_precipitation,
    BulkAlbedo,
    surface_albedo

abstract type AbstractBucketModel{FT} <: AbstractModel{FT} end
abstract type AbstractAtmosphericDrivers{FT <: AbstractFloat} end
abstract type AbstractRadiativeDrivers{FT <: AbstractFloat} end
abstract type AbstractLandAlbedoModel{FT <: AbstractFloat} end

"""
    BulkAlbedo{FT} <: AbstractLandAlbedoModel

An albedo model where the albedo of different surface types
is specified. Snow albedo is treated as constant across snow
location and across wavelength. Soil albedo is specified as a function
of latitude and longitude, but is also treated as constant across
wavelength.
"""
struct BulkAlbedo{FT} <: AbstractLandAlbedoModel{FT}
    α_snow::FT
    α_soil::Function
end



ClimaLSM.name(::AbstractBucketModel) = :bucket

"""
    struct BucketModelParameters{
        FT <: AbstractFloat,
        PSE <: AbstractEarthParameterSet,
    }

Container for holding the parameters of the bucket model.

$(DocStringExtensions.FIELDS)
"""
struct BucketModelParameters{
    FT <: AbstractFloat,
    AAM <: AbstractLandAlbedoModel,
    PSE <: AbstractEarthParameterSet,
}
    "Depth of the soil column (m) used in heat equation"
    d_soil::FT
    "Temperature at the bottom of the soil column (K); constant"
    T0::FT
    "Conductivity of the soil (W/K/m); constant"
    κ_soil::FT
    "Volumetric heat capacity of the soil (J/m^3/K); constant"
    ρc_soil::FT
    "Albedo Model"
    albedo::AAM
    "Critical SWE amount (m) where surface transitions from soil to snow"
    S_c::FT
    "Capacity of the land bucket (m)"
    W_f::FT
    "Roughness length for momentum (m)"
    z_0m::FT
    "Roughness length for scalars (m)"
    z_0b::FT
    "Earth Parameter set; physical constants, etc"
    earth_param_set::PSE
end

BucketModelParameters(
    d_soil::FT,
    T0::FT,
    κ_soil::FT,
    ρc_soil::FT,
    albedo::AAM,
    S_c::FT,
    W_f::FT,
    z_0m::FT,
    z_0b::FT,
    earth_param_set::PSE,
) where {FT, AAM, PSE} = BucketModelParameters{FT, AAM, PSE}(
    d_soil,
    T0,
    κ_soil,
    ρc_soil,
    albedo,
    S_c,
    W_f,
    z_0m,
    z_0b,
    earth_param_set,
)

"""
    PrescribedAtmosphere{FT, LP, TA, UA, QA, RA} <: AbstractAtmosphericDrivers{FT}

Container for holding prescribed atmospheric drivers and other
information needed for computing turbulent surface fluxes when
driving the bucket model in standalone mode.
$(DocStringExtensions.FIELDS)
"""
struct PrescribedAtmosphere{FT, LP, TA, UA, QA, RA} <:
       AbstractAtmosphericDrivers{FT}
    "Precipitation (m/s) function of time: positive by definition"
    liquid_precip::LP
    "Prescribed atmospheric temperature (function of time)  at the reference height (K)"
    T_atmos::TA
    "Prescribed wind speed (function of time)  at the reference height (m/s)"
    u_atmos::UA
    "Prescribed specific humidity (function of time)  at the reference height (_)"
    q_atmos::QA
    "Prescribed air density (function of time)  at the reference height (kg/m^3)"
    ρ_atmos::RA
    "Reference height, relative to surface elevation(m)"
    h_atmos::FT
    "Surface air density (kg/m^3). Note that convergence errors result if ρ_sfc = ρ_atmos."
    ρ_sfc::FT # Eventually, computed from mean surface pressure (from Atmos) + land T_sfc
end

function PrescribedAtmosphere(
    precipitation,
    T_atmos,
    u_atmos,
    q_atmos,
    ρ_atmos,
    h_atmos,
    ρ_sfc,
)
    args = (precipitation, T_atmos, u_atmos, q_atmos, ρ_atmos)
    PrescribedAtmosphere{typeof(h_atmos), typeof.(args)...}(
        args...,
        h_atmos,
        ρ_sfc,
    )
end



"""
    PrescribedRadiativeFluxes{FT, SW, LW} <: AbstractRadiativeDrivers{FT}

Container for the prescribed radiation functions needed to drive the
bucket model in standalone mode.
$(DocStringExtensions.FIELDS)
"""
struct PrescribedRadiativeFluxes{FT, SW, LW} <: AbstractRadiativeDrivers{FT}
    "Downward shortwave radiation function of time (W/m^2): positive indicates towards surface"
    SW_d::SW
    "Downward longwave radiation function of time (W/m^2): positive indicates towards surface"
    LW_d::LW
end

PrescribedRadiativeFluxes(FT, SW_d, LW_d) =
    PrescribedRadiativeFluxes{FT, typeof(SW_d), typeof(LW_d)}(SW_d, LW_d)


"""

    struct BucketModel{
         FT,
         PS <: BucketModelParameters{FT},
         ATM <: AbstractAtmosphericDrivers{FT},
         RAD <: AbstractRadiativeDrivers{FT},
         D,
     } <: AbstractBucketModel{FT}

Concrete type for the BucketModel, which store the model
domain and parameters, as well as the necessary atmosphere
and radiation fields for driving the model.
$(DocStringExtensions.FIELDS)
"""
struct BucketModel{
    FT,
    PS <: BucketModelParameters{FT},
    ATM <: AbstractAtmosphericDrivers{FT},
    RAD <: AbstractRadiativeDrivers{FT},
    D,
} <: AbstractBucketModel{FT}
    "Parameters required by the bucket model"
    parameters::PS
    "The atmospheric drivers: Prescribed or Coupled"
    atmos::ATM
    "The radiation drivers: Prescribed or Coupled"
    radiation::RAD
    "The domain of the model"
    domain::D
end

function BucketModel(;
    parameters::BucketModelParameters{FT, PSE},
    domain::ClimaLSM.Domains.AbstractDomain,
    atmosphere::ATM,
    radiation::RAD,
) where {FT, PSE, ATM, RAD}
    args = (parameters, atmosphere, radiation, domain)
    BucketModel{FT, typeof.(args)...}(args...)
end

prognostic_types(::BucketModel{FT}) where {FT} = (FT, FT, FT, FT)
prognostic_vars(::BucketModel) = (:W, :T_sfc, :Ws, :S)
auxiliary_types(::BucketModel{FT}) where {FT} = (FT, FT, FT, FT, FT)
auxiliary_vars(::BucketModel) = (:q_sfc, :E, :LHF, :SHF, :R_n)

"""
    surface_fluxes(Y,p,
                    t::FT,
                    parameters::P,
                    atmos::PA,
                    radiation::PR,
                    ) where {FT <: AbstractFloat, P <: BucketModelParameters{FT},  PA <: PrescribedAtmosphere{FT}, PR <: PrescribedRadiativeFluxes{FT}}

Computes the surface flux terms at the ground for a standalone simulation:
net radiation,  SHF,  LHF,
as well as the water vapor flux (in units of m^3/m^2/s of water).
Positive fluxes indicate flow from the ground to the atmosphere.

It solves for these given atmospheric conditions, stored in `atmos`,
 downwelling shortwave and longwave radiation (in `radiation`),
model parameters, and the surface temperature and surface specific
humidity.

Currently, we only support soil covered surfaces.
"""
function surface_fluxes(
    Y,
    p,
    t::FT,
    parameters::P,
    atmos::PA,
    radiation::PR,
) where {
    FT <: AbstractFloat,
    P <: BucketModelParameters{FT},
    PA <: PrescribedAtmosphere{FT},
    PR <: PrescribedRadiativeFluxes{FT},
}

    return surface_fluxes_at_a_point.(
        Y.bucket.T_sfc,
        p.bucket.q_sfc,
        Y.bucket.S,
        coordinate_field(Y.bucket.S),
        t,
        Ref(parameters),
        Ref(atmos),
        Ref(radiation),
    )
end

function surface_fluxes_at_a_point(
    T_sfc::FT,
    q_sfc::FT,
    S::FT,
    coords,
    t::FT,
    parameters::P,
    atmos::PA,
    radiation::PR,
) where {
    FT <: AbstractFloat,
    P <: BucketModelParameters{FT},
    PA <: PrescribedAtmosphere{FT},
    PR <: PrescribedRadiativeFluxes{FT},
}
    @unpack ρ_atmos, T_atmos, u_atmos, q_atmos, h_atmos, ρ_sfc = atmos
    @unpack albedo, z_0m, z_0b, S_c, earth_param_set = parameters
    @unpack LW_d, SW_d = radiation
    _σ = Stefan()
    _ρ_liq = ρ_cloud_liq(earth_param_set)

    # call surface fluxes for E, SHF, LHF
    ts_sfc = Thermodynamics.PhaseEquil_ρTq(earth_param_set, ρ_sfc, T_sfc, q_sfc)
    ts_in = Thermodynamics.PhaseEquil_ρTq(
        earth_param_set,
        ρ_atmos(t),
        T_atmos(t),
        q_atmos(t),
    )

    # h_atmos is relative to surface height, so we can set surface height to zero.
    state_sfc = SurfaceFluxes.SurfaceValues(FT(0), SVector{2, FT}(0, 0), ts_sfc)
    state_in = SurfaceFluxes.InteriorValues(
        h_atmos,
        SVector{2, FT}(u_atmos(t), 0),
        ts_in,
    )

    # State containers
    sc = SurfaceFluxes.ValuesOnly{FT}(;
        state_in,
        state_sfc,
        z0m = z_0m,
        z0b = z_0b,
    )

    conditions = SurfaceFluxes.surface_conditions(
        earth_param_set,
        sc,
        SurfaceFluxes.UniversalFunctions.Businger(),
    )

    α = surface_albedo(albedo, coords, S, S_c)
    # Recall that the user passed the LW and SW downwelling radiation,
    # where positive values indicate toward surface, so we need a negative sign out front
    # in order to inidicate positive R_n  = towards atmos.
    R_n = -((FT(1) - α) * SW_d(t) + LW_d(t) - _σ * T_sfc^FT(4.0))
    # Land needs a volume flux of water, not mass flux
    E = SurfaceFluxes.evaporation(sc, earth_param_set, conditions.Ch) / _ρ_liq
    return (R_n = R_n, LHF = conditions.lhf, SHF = conditions.shf, E = E)
end



"""
    make_rhs(model::BucketModel{FT}) where {FT}

Creates the rhs! function for the bucket model.
"""
function make_rhs(model::BucketModel{FT}) where {FT}
    function rhs!(dY, Y, p, t)
        @unpack d_soil, T0, κ_soil, ρc_soil, S_c, W_f = model.parameters
        # Always positive
        liquid_precip = liquid_precipitation(p, model.atmos, t)
        @unpack SHF, LHF, R_n, E = p.bucket
        F_surf_soil = @. (R_n + LHF + SHF) # Eqn 16. TODO: modify for snow

        F_bot_soil = @. -κ_soil / (d_soil) * (Y.bucket.T_sfc - T0) # Equation (16)

        E_surf_soil = @. (FT(1.0) - heaviside(Y.bucket.S)) * E # Equation (11) assuming E is volume flux
        space = axes(Y.bucket.S)

        # Always positive
        snow_melt =
            zeros(axes(Y.bucket.S)) .* (FT(1.0) .- heaviside.(Y.bucket.S)) # Equation (8)

        infiltration =
            infiltration_at_point.(
                Y.bucket.W,
                snow_melt,
                liquid_precip,
                E_surf_soil,
                W_f,
            )
        # Positive infiltration -> net (negative) flux into soil
        dY.bucket.W .= infiltration # Equation (4a) of the text. 
        dY.bucket.Ws .=
            (liquid_precip .+ snow_melt .- E_surf_soil) .- infiltration # Equation (5) of the text
        dY.bucket.S .= zeros(space) # To be equation (6)

        @. dY.bucket.T_sfc =
            -FT(1.0) / ρc_soil / d_soil * (F_surf_soil - F_bot_soil) # Eq 16
    end
    return rhs!
end
function liquid_precipitation(p, atmos::PrescribedAtmosphere, t)
    return atmos.liquid_precip(t)
end

function surface_air_density(p, atmos::PrescribedAtmosphere)
    return atmos.ρ_sfc
end

"""
    make_update_aux(model::BucketModel{FT}) where {FT}

Creates the update_aux! function for the BucketModel.
"""
function make_update_aux(model::BucketModel{FT}) where {FT}
    function update_aux!(p, Y, t)
        ρ_sfc = surface_air_density(p, model.atmos)
        p.bucket.q_sfc .=
            β.(Y.bucket.W, model.parameters.W_f) .*
            q_sat.(Y.bucket.T_sfc, Y.bucket.S, ρ_sfc, Ref(model.parameters))

        fluxes = surface_fluxes(
            Y,
            p,
            t,
            model.parameters,
            model.atmos,
            model.radiation,
        )

        @. p.bucket.SHF = fluxes.SHF
        @. p.bucket.R_n = fluxes.R_n
        @. p.bucket.E = fluxes.E
    end
    return update_aux!
end
include("./bucket_parameterizations.jl")


end

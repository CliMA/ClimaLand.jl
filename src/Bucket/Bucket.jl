module Bucket
using UnPack
using DocStringExtensions
using SurfaceFluxes
using SurfaceFluxes.UniversalFunctions
using Thermodynamics

using ClimaCore
using ClimaCore.Fields: coordinate_field, level, FieldVector
using ClimaCore.Operators: InterpolateC2F, DivergenceF2C, GradientC2F, SetValue
using ClimaCore.Geometry: WVector
using ClimaComms
using StaticArrays: SVector

using ClimaLSM
using ClimaLSM.Regridder: MapInfo, regrid_netcdf_to_field
import ..Parameters as LSMP
import ClimaLSM.Domains: coordinates, LSMSphericalShellDomain
import ClimaLSM:
    AbstractModel,
    make_update_aux,
    make_rhs,
    prognostic_vars,
    auxiliary_vars,
    name,
    prognostic_types,
    auxiliary_types,
    initialize_vars,
    initialize,
    initialize_auxiliary,
    make_set_initial_aux_state
export BucketModelParameters,
    PrescribedAtmosphere,
    PrescribedRadiativeFluxes,
    BucketModel,
    AbstractAtmosphericDrivers,
    AbstractRadiativeDrivers,
    surface_fluxes,
    net_radiation,
    surface_air_density,
    liquid_precipitation,
    snow_precipitation,
    BulkAlbedoFunction,
    BulkAlbedoMap,
    surface_albedo,
    partition_surface_fluxes

include(joinpath(pkgdir(ClimaLSM), "src/Bucket/artifacts/artifacts.jl"))

abstract type AbstractBucketModel{FT} <: AbstractModel{FT} end

"""
     AbstractAtmosphericDrivers{FT <: AbstractFloat}

An abstract type of atmospheric drivers of the bucket model.
"""
abstract type AbstractAtmosphericDrivers{FT <: AbstractFloat} end

"""
     AbstractRadiativeDrivers{FT <: AbstractFloat}

An abstract type of radiative drivers of the bucket model.
"""
abstract type AbstractRadiativeDrivers{FT <: AbstractFloat} end
abstract type AbstractLandAlbedoModel{FT <: AbstractFloat} end

"""
    BulkAlbedoFunction{FT} <: AbstractLandAlbedoModel

An albedo model where the albedo of different surface types
is specified. Snow albedo is treated as constant across snow
location and across wavelength. Surface albedo (sfc) 
is specified as a function
of latitude and longitude, but is also treated as constant across
wavelength; surface is this context refers to soil and vegetation.
"""
struct BulkAlbedoFunction{FT} <: AbstractLandAlbedoModel{FT}
    α_snow::FT
    α_sfc::Function
end

"""
    BulkAlbedoMap{FT} <: AbstractLandAlbedoModel

An albedo model where the albedo of different surface types
is specified. Snow albedo is treated as constant across snow
location and across wavelength. Surface albedo is specified via a
NetCDF file, which can be a function of time, but is treated
as constant across wavelengths; surface is this context refers 
to soil and vegetation.

Note that this option should only be used with global simulations,
i.e. with a `ClimaLSM.LSMSphericalShellDomain.`
"""
struct BulkAlbedoMap{FT} <: AbstractLandAlbedoModel{FT}
    α_snow::FT
    α_sfc::MapInfo
end

"""
    BulkAlbedoMap{FT}(regrid_dirpath;α_snow = FT(0.8), comms = ClimaComms.SingletonCommsContext()) where {FT}

Constructor for the BulkAlbedoMap that implements a default albedo map, `comms` context, and value for `α_snow`.
The `varname` must correspond to the name of the variable in the NetCDF file specified by `path`.
"""
function BulkAlbedoMap{FT}(
    regrid_dirpath;
    α_snow = FT(0.8),
    comms = ClimaComms.SingletonCommsContext(),
    path = bareground_albedo_dataset_path(),
    varname = "sw_alb",
) where {FT}
    α_sfc = MapInfo(path, varname, regrid_dirpath, comms)
    return BulkAlbedoMap{FT}(α_snow, α_sfc)
end




ClimaLSM.name(::AbstractBucketModel) = :bucket

"""
    struct BucketModelParameters{
        FT <: AbstractFloat,
        PSE,
    }

Container for holding the parameters of the bucket model.

$(DocStringExtensions.FIELDS)
"""
struct BucketModelParameters{
    FT <: AbstractFloat,
    AAM <: AbstractLandAlbedoModel,
    PSE,
}
    "Conductivity of the soil (W/K/m); constant"
    κ_soil::FT
    "Volumetric heat capacity of the soil (J/m^3/K); constant"
    ρc_soil::FT
    "Albedo Model"
    albedo::AAM
    "Critical σSWE amount (m) where surface transitions from to snow-covered"
    σS_c::FT
    "Capacity of the land bucket (m)"
    W_f::FT
    "Roughness length for momentum (m)"
    z_0m::FT
    "Roughness length for scalars (m)"
    z_0b::FT
    "τc timescale on which snow melts"
    τc::FT
    "Earth Parameter set; physical constants, etc"
    earth_param_set::PSE
end

BucketModelParameters(
    κ_soil::FT,
    ρc_soil::FT,
    albedo::AAM,
    σS_c::FT,
    W_f::FT,
    z_0m::FT,
    z_0b::FT,
    τc::FT,
    earth_param_set::PSE,
) where {FT, AAM, PSE} = BucketModelParameters{FT, AAM, PSE}(
    κ_soil,
    ρc_soil,
    albedo,
    σS_c,
    W_f,
    z_0m,
    z_0b,
    τc,
    earth_param_set,
)

"""
    PrescribedAtmosphere{FT, LP, SP, TA, UA, QA, RA} <: AbstractAtmosphericDrivers{FT}

Container for holding prescribed atmospheric drivers and other
information needed for computing turbulent surface fluxes when
driving the bucket model in standalone mode.
$(DocStringExtensions.FIELDS)
"""
struct PrescribedAtmosphere{FT, LP, SP, TA, UA, QA, RA} <:
       AbstractAtmosphericDrivers{FT}
    "Precipitation (m/s) function of time: positive by definition"
    liquid_precip::LP
    "Snow precipitation (m/s) function of time: positive by definition"
    snow_precip::SP
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
    snow_precipitation,
    T_atmos,
    u_atmos,
    q_atmos,
    ρ_atmos,
    h_atmos,
    ρ_sfc,
)
    args =
        (precipitation, snow_precipitation, T_atmos, u_atmos, q_atmos, ρ_atmos)
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
    domain::ClimaLSM.Domains.AbstractLSMDomain,
    atmosphere::ATM,
    radiation::RAD,
) where {FT, PSE, ATM, RAD}
    if parameters.albedo isa BulkAlbedoMap
        typeof(domain) <: LSMSphericalShellDomain ? nothing :
        error("Using an albedo map requires a global run.")
    end

    args = (parameters, atmosphere, radiation, domain)
    BucketModel{FT, typeof.(args)...}(args...)
end



prognostic_types(::BucketModel{FT}) where {FT} = (FT, FT, FT, FT)
prognostic_vars(::BucketModel) = (:W, :T, :Ws, :σS)
auxiliary_types(::BucketModel{FT}) where {FT} = (FT, FT, FT, FT, FT, FT)
auxiliary_vars(::BucketModel) =
    (:q_sfc, :evaporation, :turbulent_energy_flux, :R_n, :T_sfc, :α_sfc)

"""
    ClimaLSM.initialize(model::BucketModel{FT}) where {FT}

Initializes the variables for the `BucketModel`.

Note that the `BucketModel` has prognostic variables that are defined on different
subsets of the domain. Because of that, we have to treat them independently.
In LSM models which are combinations of standalone component models, this is not 
needed, and we can use the default `initialize`. Here, however, we need to do some
hardcoding specific to this model.
"""
function ClimaLSM.initialize(model::BucketModel{FT}) where {FT}
    model_name = name(model)
    subsurface_coords, surface_coords =
        ClimaLSM.Domains.coordinates(model.domain)
    # Temperature `T` is the only prognostic variable on the subsurface.
    subsurface_prog =
        ClimaLSM.initialize_vars((:T,), (FT,), subsurface_coords, model_name)

    # Surface variables:
    surface_keys = [key for key in prognostic_vars(model) if key != :T]
    surface_types = [FT for _ in surface_keys]
    surface_prog = ClimaLSM.initialize_vars(
        surface_keys,
        surface_types,
        surface_coords,
        model_name,
    )
    surface_prog_states = map(surface_keys) do (key)
        getproperty(surface_prog.bucket, key)
    end

    values = (surface_prog_states..., subsurface_prog.bucket.T)
    keys = (surface_keys..., :T)

    Y = ClimaCore.Fields.FieldVector(; model_name => (; zip(keys, values)...))

    # All aux variables for this model live on the surface
    p = initialize_auxiliary(model, surface_coords)

    return Y, p, ClimaLSM.Domains.coordinates(model.domain)
end

"""
    ClimaLSM.make_set_initial_aux_state(model::BucketModel{FT}) where{FT}

Returns the set_initial_aux_state! function, which updates the auxiliary
state `p` in place with the initial values corresponding to Y(t=t0) = Y0.

In this case, we also use this method to update the initial values for the
spatially varying parameter fields, read in from data files.
"""
function ClimaLSM.make_set_initial_aux_state(model::BucketModel)
    update_aux! = make_update_aux(model)
    function set_initial_aux_state!(p, Y0, t0)
        set_initial_parameter_field!(
            model.parameters.albedo,
            p,
            ClimaCore.Fields.coordinate_field(model.domain.surface.space),
        )
        update_aux!(p, Y0, t0)
    end
    return set_initial_aux_state!
end

"""
    function set_initial_parameter_field!(
        albedo::BulkAlbedoFunction{FT},
        p,
        surface_coords,
    ) where {FT}

Updates the spatially-varying but constant in time surface
 albedo stored in the
auxiliary vector `p` in place,  according to the
passed function of latitute and longitude stored in `albedo.α_sfc`.
"""
function set_initial_parameter_field!(
    albedo::BulkAlbedoFunction{FT},
    p,
    surface_coords,
) where {FT}
    (; α_sfc) = albedo
    @. p.bucket.α_sfc = α_sfc(surface_coords)
end

"""
    function set_initial_parameter_field!(
        albedo::BulkAlbedoMap{FT},
        p,
        surface_coords,
    ) where {FT}

Updates spatially-varying surface albedo stored in the
auxiliary vector `p` in place, according to a
NetCDF file. This data file is encapsulated in an object of
type `ClimaLSM.Regridder.MapInfo` in the field albedo.α_sfc.

The NetCDF file is read in, regridded, and projected onto 
the surface space of the LSM using ClimaCoreTempestRemap. The result
is a ClimaCore.Fields.Field of albedo values.
"""
function set_initial_parameter_field!(
    albedo::BulkAlbedoMap{FT},
    p,
    surface_coords,
) where {FT}
    α_sfc = albedo.α_sfc
    (; comms, varname, path, regrid_dirpath) = α_sfc
    p.bucket.α_sfc .= regrid_netcdf_to_field(
        FT,
        regrid_dirpath,
        comms,
        path,
        varname,
        axes(surface_coords),
    )
end



"""
    surface_fluxes(Y,p,
                    t::FT,
                    parameters::P,
                    atmos::PA,
                    ) where {FT <: AbstractFloat, P <: BucketModelParameters{FT},  PA <: PrescribedAtmosphere{FT}}

Computes the turbulent surface flux terms at the ground for a standalone simulation:
turbulent energy fluxes
as well as the water vapor flux (in units of m^3/m^2/s of water).
Positive fluxes indicate flow from the ground to the atmosphere.

It solves for these given atmospheric conditions, stored in `atmos`,
model parameters, and the surface temperature and surface specific
humidity.

"""
function surface_fluxes(
    Y,
    p,
    t::FT,
    parameters::P,
    atmos::PA,
) where {
    FT <: AbstractFloat,
    P <: BucketModelParameters{FT},
    PA <: PrescribedAtmosphere{FT},
}
    β_sfc = beta_factor.(Y.bucket.W, Y.bucket.σS, parameters.W_f)
    return surface_fluxes_at_a_point.(
        p.bucket.T_sfc,
        p.bucket.q_sfc,
        t,
        β_sfc,
        Ref(parameters),
        Ref(atmos),
    )
end

function surface_fluxes_at_a_point(
    T_sfc::FT,
    q_sfc::FT,
    t::FT,
    β_sfc::FT,
    parameters::P,
    atmos::PA,
) where {
    FT <: AbstractFloat,
    P <: BucketModelParameters{FT},
    PA <: PrescribedAtmosphere{FT},
}
    @unpack ρ_atmos, T_atmos, u_atmos, q_atmos, h_atmos, ρ_sfc = atmos
    @unpack z_0m, z_0b, earth_param_set = parameters
    _ρ_liq = LSMP.ρ_cloud_liq(earth_param_set)

    thermo_params = LSMP.thermodynamic_parameters(earth_param_set)

    # call surface fluxes for E, energy fluxes
    ts_sfc = Thermodynamics.PhaseEquil_ρTq(thermo_params, ρ_sfc, T_sfc, q_sfc)
    ts_in = Thermodynamics.PhaseEquil_ρTq(
        thermo_params,
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
        beta = β_sfc,
    )
    surface_flux_params = LSMP.surface_fluxes_parameters(earth_param_set)
    conditions = SurfaceFluxes.surface_conditions(surface_flux_params, sc)

    # Land needs a volume flux of water, not mass flux
    evaporation =
        SurfaceFluxes.evaporation(surface_flux_params, sc, conditions.Ch) /
        _ρ_liq
    return (
        turbulent_energy_flux = conditions.lhf .+ conditions.shf,
        evaporation = evaporation,
    )
end



"""
    net_radiation(
        Y,
        p,
        t::FT,
        parameters::P,
        radiation::PR,
        ) where {FT <: AbstractFloat, P <: BucketModelParameters{FT},  PA <: PrescribedAtmosphere{FT}}

Computes the net radiation at the ground for a standalone simulation.
Positive fluxes indicate flow from the ground to the atmosphere.
"""
function net_radiation(
    Y,
    p,
    t::FT,
    parameters::P,
    radiation::PR,
) where {
    FT <: AbstractFloat,
    P <: BucketModelParameters{FT},
    PR <: PrescribedRadiativeFluxes{FT},
}
    return radiative_fluxes_at_a_point.(
        p.bucket.T_sfc,
        p.bucket.α_sfc,
        Y.bucket.σS,
        t,
        Ref(parameters),
        Ref(radiation),
    )
end

function radiative_fluxes_at_a_point(
    T_sfc::FT,
    α_sfc::FT,
    σS::FT,
    t::FT,
    parameters::P,
    radiation::PR,
) where {
    FT <: AbstractFloat,
    P <: BucketModelParameters{FT},
    PR <: PrescribedRadiativeFluxes{FT},
}
    @unpack albedo, earth_param_set, σS_c = parameters
    @unpack LW_d, SW_d = radiation
    _σ = LSMP.Stefan(earth_param_set)
    (; α_snow) = albedo
    α = surface_albedo(α_sfc, α_snow, σS, σS_c)
    # Recall that the user passed the LW and SW downwelling radiation,
    # where positive values indicate toward surface, so we need a negative sign out front
    # in order to inidicate positive R_n  = towards atmos.
    R_n = -(1 - α) * SW_d(t) - LW_d(t) + _σ * T_sfc^4
    return R_n
end


"""
    make_rhs(model::BucketModel{FT}) where {FT}

Creates the rhs! function for the bucket model.
"""
function make_rhs(model::BucketModel{FT}) where {FT}
    function rhs!(dY, Y, p, t)
        @unpack κ_soil, ρc_soil, σS_c, W_f = model.parameters

        #Currently, the entire surface is assumed to be
        # snow covered entirely or not at all.
        snow_cover_fraction = heaviside.(Y.bucket.σS)

        # In this case, there is just one set of surface fluxes to compute.
        # Since q_sfc itself has already been modified to account for
        # snow covered regions, and since the surface temperature is
        # assumed to be the same for snow and underlying land,
        # the returned fluxes are computed correctly for the cell
        # regardless of snow-coverage.

        # The below is NOT CORRECT if we want the snow
        # cover fraction to be intermediate between 0 and 1.
        @unpack turbulent_energy_flux, R_n, evaporation = p.bucket
        F_sfc = @. (R_n + turbulent_energy_flux) # Eqn (15)

        _T_freeze = LSMP.T_freeze(model.parameters.earth_param_set)
        _LH_f0 = LSMP.LH_f0(model.parameters.earth_param_set)
        _ρ_liq = LSMP.ρ_cloud_liq(model.parameters.earth_param_set)
        _ρLH_f0 = _ρ_liq * _LH_f0 # Latent heat per unti volume.
        # partition energy fluxes
        partitioned_fluxes =
            partition_surface_fluxes.(
                Y.bucket.σS,
                p.bucket.T_sfc,
                model.parameters.τc,
                snow_cover_fraction,
                evaporation,
                F_sfc,
                _ρLH_f0,
                _T_freeze,
            )
        (; F_melt, F_into_snow, G) = partitioned_fluxes

        # Temperature profile of soil. 
        gradc2f = ClimaCore.Operators.GradientC2F()
        divf2c = ClimaCore.Operators.DivergenceF2C(
            top = ClimaCore.Operators.SetValue(ClimaCore.Geometry.WVector.(G)),
            bottom = ClimaCore.Operators.SetValue(
                ClimaCore.Geometry.WVector.(FT(0.0)),
            ),
        )
        @. dY.bucket.T = -1 / ρc_soil * (divf2c(-κ_soil * gradc2f(Y.bucket.T))) # Simple heat equation, Eq 6


        # Partition water fluxes
        liquid_precip = liquid_precipitation(p, model.atmos, t) # always positive
        snow_precip = snow_precipitation(p, model.atmos, t) # always positive
        # Always positive; F_melt at present already has σ factor in it.
        snow_melt = @. (-F_melt / _ρLH_f0) # Equation (20)

        infiltration = @. infiltration_at_point(
            Y.bucket.W,
            snow_melt, # snow melt already multipied by snow_cover_fraction
            liquid_precip,
            (1 - snow_cover_fraction) * evaporation,
            W_f,
        ) # Equation (4b) of the text.

        # Positive infiltration -> net (negative) flux into soil
        dY.bucket.W .= infiltration # Equation (4) of the text.

        dY.bucket.Ws = @. (
            (
                liquid_precip + snow_melt -
                (1 - snow_cover_fraction) * evaporation
            ) - infiltration
        ) # Equation (5) of the text, snow melt already multipied by snow_cover_fraction

        # snow melt already multipied by snow_cover_fraction
        dY.bucket.σS =
            @. (snow_precip - snow_cover_fraction * evaporation - snow_melt) # Equation (11)
    end
    return rhs!
end

"""
    liquid_precipitation(p, atmos::PrescribedAtmosphere, t)

Returns the liquid precipitation (m/s) at the surface.
"""
function liquid_precipitation(p, atmos::PrescribedAtmosphere, t)
    return atmos.liquid_precip(t)
end

"""
    snow_precipitation(p, atmos::PrescribedAtmosphere, t)

Returns the precipitation in snow (m of liquid water/s) at the surface.
"""
function snow_precipitation(p, atmos::PrescribedAtmosphere, t)
    return atmos.snow_precip(t)
end

"""
    surface_air_density(p, atmos::PrescribedAtmosphere)

Returns the air density (kg/m^3) at the surface.
"""
function surface_air_density(p, atmos::PrescribedAtmosphere)
    return atmos.ρ_sfc
end

"""
    make_update_aux(model::BucketModel{FT}) where {FT}

Creates the update_aux! function for the BucketModel.
"""
function make_update_aux(model::BucketModel{FT}) where {FT}
    function update_aux!(p, Y, t)
        face_space =
            ClimaLSM.Domains.obtain_face_space(model.domain.subsurface.space)
        N = ClimaCore.Spaces.nlevels(face_space)
        interp_c2f = ClimaCore.Operators.InterpolateC2F(
            top = ClimaCore.Operators.Extrapolate(),
            bottom = ClimaCore.Operators.Extrapolate(),
        )
        surface_space = model.domain.surface.space
        p.bucket.T_sfc .= ClimaCore.Fields.Field(
            ClimaCore.Fields.field_values(
                ClimaCore.Fields.level(
                    interp_c2f.(Y.bucket.T),
                    ClimaCore.Utilities.PlusHalf(N - 1),
                ),
            ),
            surface_space,
        )
        ρ_sfc = surface_air_density(p, model.atmos)

        # This relies on the surface specific humidity being computed
        # entirely over snow or over soil. i.e. the snow cover fraction must be a heaviside
        # here, otherwise we would need two values of q_sfc!
        p.bucket.q_sfc .=
            saturation_specific_humidity.(
                p.bucket.T_sfc,
                Y.bucket.σS,
                ρ_sfc,
                Ref(model.parameters),
            )

        turbulent_fluxes =
            surface_fluxes(Y, p, t, model.parameters, model.atmos)

        p.bucket.turbulent_energy_flux .= turbulent_fluxes.turbulent_energy_flux
        p.bucket.evaporation .= turbulent_fluxes.evaporation
        p.bucket.R_n .=
            net_radiation(Y, p, t, model.parameters, model.radiation)

    end
    return update_aux!
end
include("./bucket_parameterizations.jl")


end

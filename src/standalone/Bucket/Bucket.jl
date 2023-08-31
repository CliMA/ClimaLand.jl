module Bucket
using UnPack
using DocStringExtensions
using Thermodynamics
using Dates
using NCDatasets
using ClimaCore
using ClimaCore.Fields: coordinate_field, level, FieldVector
using ClimaCore.Operators: InterpolateC2F, DivergenceF2C, GradientC2F, SetValue
using ClimaCore.Geometry: WVector
using ClimaComms

using ClimaLSM

using ClimaLSM.Regridder: MapInfo, regrid_netcdf_to_field
using ClimaLSM.FileReader
import ..Parameters as LSMP
import ClimaLSM.Domains: coordinates, LSMSphericalShellDomain
using ClimaLSM:
    AbstractAtmosphericDrivers,
    AbstractRadiativeDrivers,
    liquid_precipitation,
    snow_precipitation,
    surface_fluxes,
    net_radiation,
    construct_atmos_ts,
    compute_ρ_sfc,
    AbstractExpModel,
    heaviside,
    PrescribedAtmosphere,
    add_dss_buffer_to_aux
import ClimaLSM:
    make_update_aux,
    make_compute_exp_tendency,
    prognostic_vars,
    auxiliary_vars,
    name,
    prognostic_types,
    auxiliary_types,
    initialize_vars,
    initialize,
    initialize_auxiliary,
    make_set_initial_aux_state,
    surface_temperature,
    surface_air_density,
    surface_specific_humidity,
    surface_evaporative_scaling,
    surface_albedo,
    surface_emissivity,
    surface_height
export BucketModelParameters,
    BucketModel,
    BulkAlbedoFunction,
    BulkAlbedoStatic,
    BulkAlbedoTemporal,
    surface_albedo,
    partition_surface_fluxes

include(
    joinpath(pkgdir(ClimaLSM), "src/standalone/Bucket/artifacts/artifacts.jl"),
)

abstract type AbstractBucketModel{FT} <: AbstractExpModel{FT} end

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
    BulkAlbedoStatic{FT} <: AbstractLandAlbedoModel

An albedo model where the albedo of different surface types
is specified. Snow albedo is treated as constant across snow
location and across wavelength. Surface albedo is specified via a
NetCDF file, which can be a function of time, but is treated
as constant across wavelengths; surface is this context refers
to soil and vegetation. This albedo type is static in time.

Note that this option should only be used with global simulations,
i.e. with a `ClimaLSM.LSMSphericalShellDomain.`
"""
struct BulkAlbedoStatic{FT} <: AbstractLandAlbedoModel{FT}
    α_snow::FT
    α_sfc::MapInfo
end

"""
    BulkAlbedoStatic{FT}(
        regrid_dirpath::String;
        α_snow = FT(0.8),
        comms = ClimaComms.SingletonCommsContext(),
        varname = "sw_alb",
        path = Bucket.bareground_albedo_dataset_path(),
    ) where {FT}

Constructor for the BulkAlbedoStatic that implements a default albedo map, `comms` context, and value for `α_snow`.
The `varname` must correspond to the name of the variable in the NetCDF file specified by `path`.

The `bareground_albedo_dataset_path` artifact can be used as a default with this type.
"""
function BulkAlbedoStatic{FT}(
    regrid_dirpath::String;
    α_snow = FT(0.8),
    comms = ClimaComms.SingletonCommsContext(),
    varname = "sw_alb",
    path = Bucket.bareground_albedo_dataset_path(),
) where {FT}
    α_sfc = MapInfo(path, varname, regrid_dirpath, comms)
    return BulkAlbedoStatic{FT}(α_snow, α_sfc)
end

"""
    BulkAlbedoTemporal{FT} <: AbstractLandAlbedoModel

An albedo model where the albedo of different surface types
is specified. Albedo is specified via a NetCDF file which is a function
of time and covers all surface types (soil, vegetation, snow, etc).
This albedo type changes over time according to the input file.

Note that this option should only be used with global simulations,
i.e. with a `ClimaLSM.LSMSphericalShellDomain.`
"""
struct BulkAlbedoTemporal{FT} <: AbstractLandAlbedoModel{FT}
    albedo_info::FileReader.PrescribedData
end

"""
    BulkAlbedoTemporal{FT}(
        regrid_dirpath::String,
        date_ref::Union{DateTime, DateTimeNoLeap},
        t_start::FT,
        Space::ClimaCore.Spaces.AbstractSpace;
        input_file = Bucket.cesm2_albedo_dataset_path(),
        varname = "sw_alb"
    ) where {FT}

Constructor for the BulkAlbedoTemporal struct.
The `varname` must correspond to the name of the variable in the NetCDF
file specified by `input_file`.
The input data file must have a time component; otherwise BulkAlbedoStatic
should be used.
"""
function BulkAlbedoTemporal{FT}(
    regrid_dirpath::String,
    date_ref::Union{DateTime, DateTimeNoLeap},
    t_start::FT,
    space::ClimaCore.Spaces.AbstractSpace;
    input_file = Bucket.cesm2_albedo_dataset_path(),
    varname = "sw_alb",
) where {FT}
    # Verify inputs
    if typeof(space) <: ClimaCore.Spaces.PointSpace
        error("Using an albedo map requires a global run.")
    end
    NCDataset(input_file, "r") do ds
        if !("time" in keys(ds))
            error(
                "Using a temporal albedo map requires data with time dimension.",
            )
        end
    end

    # Construct object containing info to read in surface albedo over time
    data_info = prescribed_data_init(
        regrid_dirpath,
        input_file,
        varname,
        date_ref,
        t_start,
        space,
    )
    return BulkAlbedoTemporal{FT}(data_info)
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

"""
   BucketModel(; parameters::BucketModelParameters{FT, PSE},
                 domain::D,
                 atmosphere::ATM,
                 radiation::RAD,
               ) where {FT, PSE, ATM, RAD, D<: ClimaLSM.Domains.AbstractLSMDomain}

An outer constructor for the `BucketModel`, which enforces the
constraints:
1. The bucket model domain is of type <: ClimaLSM.Domains.AbstractLSMDomain
2. Using an albedo read from a lat/lon file requires a global run.

Since the bucket model has prognostic variables that live on the surface of the domain
(snow water equivalent, bucket water content, runoff), as well as defined within
the entire domain (temperature as a function of depth), we make use of
ClimaLSM.Domains.AbstractLSMDomains which are set up for this use case.
"""
function BucketModel(;
    parameters::BucketModelParameters{FT, PSE},
    domain::D,
    atmosphere::ATM,
    radiation::RAD,
) where {FT, PSE, ATM, RAD, D <: ClimaLSM.Domains.AbstractLSMDomain}
    if parameters.albedo isa Union{BulkAlbedoStatic, BulkAlbedoTemporal}
        typeof(domain) <: LSMSphericalShellDomain ? nothing :
        error("Using an albedo map requires a global run.")
    end

    args = (parameters, atmosphere, radiation, domain)
    BucketModel{FT, typeof.(args)...}(args...)
end



prognostic_types(::BucketModel{FT}) where {FT} = (FT, FT, FT, FT)
prognostic_vars(::BucketModel) = (:W, :T, :Ws, :σS)
auxiliary_types(::BucketModel{FT}) where {FT} = (FT, FT, FT, FT, FT, FT, FT)
auxiliary_vars(::BucketModel) =
    (:q_sfc, :evaporation, :turbulent_energy_flux, :R_n, :T_sfc, :α_sfc, :ρ_sfc)

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

    p = initialize_auxiliary(model, surface_coords)
    p = add_dss_buffer_to_aux(p, model.domain.surface)
    p = add_dss_buffer_to_aux(p, model.domain.subsurface)
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
        albedo::BulkAlbedoStatic{FT},
        p,
        surface_coords,
    ) where {FT}

Initializes spatially-varying surface albedo stored in the
auxiliary vector `p` in place, according to a
NetCDF file. This data file is encapsulated in an object of
type `ClimaLSM.Regridder.MapInfo` in the field albedo.α_sfc.

The NetCDF file is read in, regridded, and projected onto
the surface space of the LSM using ClimaCoreTempestRemap. The result
is a ClimaCore.Fields.Field of albedo values.
"""
function set_initial_parameter_field!(
    albedo::BulkAlbedoStatic{FT},
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
    function set_initial_parameter_field!(
        albedo::BulkAlbedoTemporal{FT},
        p,
        surface_coords,
    ) where {FT}

Initializes spatially- and temporally-varying surface albedo stored in
the auxiliary vector `p` in place, according to a
NetCDF file. This data file is encapsulated in an object of
type `ClimaLSM.FileReader.PrescribedData` in the field albedo.albedo_info.
This object contains a reference date and start time, which are used
to get the start date.

The NetCDF file is read in at the dates closest to this start date,
regridded, and projected onto
the surface space of the LSM using ClimaCoreTempestRemap. The result
is a ClimaCore.Fields.Field of albedo values.
"""
function set_initial_parameter_field!(
    albedo::BulkAlbedoTemporal{FT},
    p,
    surface_coords,
) where {FT}
    sim_info = albedo.albedo_info.sim_info
    date_start = sim_info.date_ref + Dates.Second(round(sim_info.t_start))
    space = axes(surface_coords)

    read_data_fields!(albedo.albedo_info, date_start, space)
    p.bucket.α_sfc .= interpolate_data(albedo.albedo_info, date_start, space)
end

"""
    make_compute_exp_tendency(model::BucketModel{FT}) where {FT}

Creates the compute_exp_tendency! function for the bucket model.
"""
function make_compute_exp_tendency(model::BucketModel{FT}) where {FT}
    function compute_exp_tendency!(dY, Y, p, t)
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
        liquid_precip = liquid_precipitation(model.atmos, p, t) # always positive
        snow_precip = snow_precipitation(model.atmos, p, t) # always positive
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
    return compute_exp_tendency!
end

"""
    make_update_aux(model::BucketModel{FT}) where {FT}

Creates the update_aux! function for the BucketModel.
"""
function make_update_aux(model::BucketModel{FT}) where {FT}
    function update_aux!(p, Y, t)
        p.bucket.T_sfc .= ClimaLSM.Domains.top_center_to_surface(Y.bucket.T)
        p.bucket.ρ_sfc .=
            surface_air_density(model.atmos, model, Y, p, t, p.bucket.T_sfc)

        # This relies on the surface specific humidity being computed
        # entirely over snow or over soil. i.e. the snow cover fraction must be a heaviside
        # here, otherwise we would need two values of q_sfc!
        p.bucket.q_sfc .=
            saturation_specific_humidity.(
                p.bucket.T_sfc,
                Y.bucket.σS,
                p.bucket.ρ_sfc,
                Ref(model.parameters),
            )

        # Compute turbulent surface fluxes
        conditions = surface_fluxes(model.atmos, model, Y, p, t)
        p.bucket.turbulent_energy_flux .= conditions.lhf .+ conditions.shf
        p.bucket.evaporation .= conditions.vapor_flux

        # Radiative fluxes
        p.bucket.R_n .= net_radiation(model.radiation, model, Y, p, t)

        # Surface albedo
        p.bucket.α_sfc .=
            next_albedo(model.parameters.albedo, model.parameters, Y, p, t)
    end
    return update_aux!
end

"""
    next_albedo(model_albedo::Union{BulkAlbedoFunction{FT}, BulkAlbedoStatic{FT}},
        parameters, Y, p, t)

Update the surface albedo for time `t`. These albedo model types aren't explicitly
dependent on `t`, but depend on quantities which may change over time.

The albedo is calculated by linearly interpolating between the albedo
of snow and of the surface, based on the snow water equivalent `S` relative to
the parameter `S_c`. The linear interpolation is taken from Lague et al 2019.
"""
function next_albedo(
    model_albedo::Union{BulkAlbedoFunction{FT}, BulkAlbedoStatic{FT}},
    parameters,
    Y,
    p,
    t,
) where {FT}
    (; α_snow) = model_albedo
    (; σS_c) = parameters
    α_sfc = p.bucket.α_sfc
    σS = Y.bucket.σS
    return @. ((1 - σS / (σS + σS_c)) * α_sfc + σS / (σS + σS_c) * α_snow)
end

"""
    next_albedo(model_albedo::BulkAlbedoTemporal{FT}, parameters, Y, p, t)

Update the surface albedo for time t. For a file containing albedo
information over time, this reads in the value for time t.
"""
function next_albedo(
    model_albedo::BulkAlbedoTemporal{FT},
    parameters,
    Y,
    p,
    t,
) where {FT}
    # Get the current date from `t`
    sim_info = model_albedo.albedo_info.sim_info
    sim_date = to_datetime(
        sim_info.date_ref + Second(round(sim_info.t_start)) + Second(round(t)),
    )
    # Use next date if it's closest to current time
    # This maintains `all_dates[date_idx]` <= `sim_date` < `all_dates[date_idx + 1]`
    if sim_date >= to_datetime(next_date_in_file(model_albedo.albedo_info))
        read_data_fields!(model_albedo.albedo_info, sim_date, axes(Y.bucket.W))
    end
    # Interpolate data value to current time
    return interpolate_data(
        model_albedo.albedo_info,
        sim_date,
        axes(Y.bucket.W),
    )
end

include("./bucket_parameterizations.jl")

end

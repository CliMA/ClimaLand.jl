module Bucket
using DocStringExtensions
using Thermodynamics
using Dates
using NCDatasets
using ClimaCore
using ClimaCore.Fields: coordinate_field, level, FieldVector
using ClimaCore.Operators: InterpolateC2F, DivergenceF2C, GradientC2F, SetValue
using ClimaCore.Geometry: WVector
using ClimaComms

using ClimaLand

using ClimaLand.FileReader
import ..Parameters as LP
import ClimaLand.Domains: coordinates, SphericalShell
using ClimaLand:
    AbstractAtmosphericDrivers,
    AbstractRadiativeDrivers,
    liquid_precipitation,
    snow_precipitation,
    turbulent_fluxes,
    net_radiation,
    construct_atmos_ts,
    compute_ρ_sfc,
    AbstractExpModel,
    heaviside,
    PrescribedAtmosphere,
    add_dss_buffer_to_aux
import ClimaLand:
    make_update_aux,
    make_compute_exp_tendency,
    prognostic_vars,
    auxiliary_vars,
    name,
    prognostic_types,
    auxiliary_types,
    prognostic_domain_names,
    auxiliary_domain_names,
    initialize_vars,
    initialize,
    initialize_auxiliary,
    make_set_initial_cache,
    surface_temperature,
    surface_air_density,
    surface_evaporative_scaling,
    surface_albedo,
    surface_emissivity,
    surface_height,
    get_drivers
export BucketModelParameters,
    BucketModel,
    BulkAlbedoFunction,
    BulkAlbedoStatic,
    BulkAlbedoTemporal,
    surface_albedo,
    partition_snow_surface_fluxes

include(
    joinpath(pkgdir(ClimaLand), "src/standalone/Bucket/artifacts/artifacts.jl"),
)

abstract type AbstractBucketModel{FT} <: AbstractExpModel{FT} end

abstract type AbstractLandAlbedoModel{FT <: AbstractFloat} end

"""
    BulkAlbedoFunction{FT, F <: ClimaCore.Fields.Field} <: AbstractLandAlbedoModel

An albedo model where the albedo of different surface types
is specified. Snow albedo is treated as constant across snow
location and across wavelength. Bareground surface albedo
is specified as a function
of latitude and longitude, but is also treated as constant across
wavelength; bareground surface is this context refers to soil and vegetation.
"""
struct BulkAlbedoFunction{FT, F <: ClimaCore.Fields.Field} <:
       AbstractLandAlbedoModel{FT}
    α_snow::FT
    α_bareground::F
end

function BulkAlbedoFunction{FT}(
    α_snow::FT,
    α_bareground::F,
) where {FT <: AbstractFloat, F}
    BulkAlbedoFunction(α_snow, α_bareground)
end

function BulkAlbedoFunction{FT}(
    α_snow::FT,
    α_bareground_func::Function,
    space,
) where {FT <: AbstractFloat}
    α_bareground = SpaceVaryingInput(α_bareground_func, space)
    args = (α_snow, α_bareground)
    BulkAlbedoFunction{typeof.(args)...}(args...)
end

"""
    BulkAlbedoStatic{FT, F <: ClimaCore.Fields.Field} <: AbstractLandAlbedoModel

An albedo model where the albedo of different surface types
is specified. Snow albedo is treated as constant across snow
location and across wavelength. Bareground surface albedo is specified
via a NetCDF file that is treated
as constant across wavelengths; surface is this context refers
to soil and vegetation. This albedo type is static in time.

Note that this option should only be used with global simulations,
i.e. with a `ClimaLand.LSMSphericalShellDomain.`
"""
struct BulkAlbedoStatic{FT, F <: ClimaCore.Fields.Field} <:
       AbstractLandAlbedoModel{FT}
    α_snow::FT
    α_bareground::F
end

function BulkAlbedoStatic{FT}(α_snow::FT, α_bareground::F) where {FT, F}
    return BulkAlbedoStatic(α_snow, α_bareground)
end

"""
    BulkAlbedoStatic{FT}(
        regrid_dirpath::String,
        surface_space::ClimaCore.Spaces.AbstractSpace;
        α_snow = FT(0.8),
        varname = ["sw_alb"],
        get_infile::Function = Bucket.bareground_albedo_dataset_path,
    ) where {FT}

Constructor for the BulkAlbedoStatic that implements a default albedo map,
`comms` context, and value for `α_snow`.
The `varname` must correspond to the name of the variable in the NetCDF
file retrieved by `infile_path`.
`infile_path` is a function that uses ArtifactWrappers.jl to return a path to
the data file and download the data if it doesn't already exist on the machine.

The `bareground_albedo_dataset_path` artifact will be used as a default with this type.
"""
function BulkAlbedoStatic{FT}(
    regrid_dirpath::String,
    surface_space::ClimaCore.Spaces.AbstractSpace;
    α_snow = FT(0.8),
    varnames = ["sw_alb"],
    get_infile::Function = Bucket.bareground_albedo_dataset_path,
) where {FT}
    if surface_space isa ClimaCore.Spaces.PointSpace
        error("Using an albedo map requires a global run.")
    end
    α_bareground_data = PrescribedDataStatic{FT}(
        get_infile,
        regrid_dirpath,
        varnames,
        surface_space,
    )

    # Albedo file only has one variable, so access first `varname`
    varname = varnames[1]
    α_bareground = SpaceVaryingInput(α_bareground_data, varname, surface_space)
    return BulkAlbedoStatic(α_snow, α_bareground)
end

"""
    BulkAlbedoTemporal{FT, FR <: FileReader.PrescribedDataTemporal}
                       <: AbstractLandAlbedoModel

An albedo model where the albedo of different surface types
is specified. Albedo is specified via a NetCDF file which is a function
of time and covers all surface types (soil, vegetation, snow, etc).
This albedo type changes over time according to the input file.

Note that this option should only be used with global simulations,
i.e. with a `ClimaLand.LSMSphericalShellDomain.`
"""
struct BulkAlbedoTemporal{FT, FR <: FileReader.PrescribedDataTemporal} <:
       AbstractLandAlbedoModel{FT}
    albedo_info::FR
end

"""
    BulkAlbedoTemporal{FT}(
        regrid_dirpath::String,
        date_ref::Union{DateTime, DateTimeNoLeap},
        t_start,
        Space::ClimaCore.Spaces.AbstractSpace;
        get_infile = Bucket.cesm2_albedo_dataset_path,
        varname = "sw_alb"
    ) where {FT}

Constructor for the BulkAlbedoTemporal struct.
The `varname` must correspond to the name of the variable in the NetCDF
file retrieved by the `get_infile` function.
`get_infile` uses ArtifactWrappers.jl to return a path to the data file
and download the data if it doesn't already exist on the machine.
The input data file must have a time component; otherwise BulkAlbedoStatic
should be used.
"""
function BulkAlbedoTemporal{FT}(
    regrid_dirpath::String,
    date_ref::Union{DateTime, DateTimeNoLeap},
    t_start,
    space::ClimaCore.Spaces.AbstractSpace;
    get_infile = Bucket.cesm2_albedo_dataset_path,
    varname = "sw_alb",
) where {FT}
    # Verify inputs
    if typeof(space) <: ClimaCore.Spaces.PointSpace
        error("Using an albedo map requires a global run.")
    end

    # Construct object containing info to read in surface albedo over time
    data_info = PrescribedDataTemporal{FT}(
        regrid_dirpath,
        get_infile,
        [varname],
        date_ref,
        t_start,
        space,
    )
    return BulkAlbedoTemporal{FT, typeof(data_info)}(data_info)
end



ClimaLand.name(::AbstractBucketModel) = :bucket

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
    "Fraction of critical amount of snow at which sublimation β begins to decay to zero (unitless)"
    f_snow::FT
    "Capacity of the land bucket (m)"
    W_f::FT
    "Fraction of bucket capacity at which evaporation β begins to decay to zero (unitless)"
    f_bucket::FT
    "Exponent used in β decay (unitless)"
    p::FT
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
    earth_param_set::PSE;
    f_snow = FT(0.0),
    f_bucket = FT(0.75),
    p = FT(1),
) where {FT, AAM, PSE} = BucketModelParameters{FT, AAM, PSE}(
    κ_soil,
    ρc_soil,
    albedo,
    σS_c,
    f_snow,
    W_f,
    f_bucket,
    p,
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
               ) where {FT, PSE, ATM, RAD, D<: ClimaLand.Domains.AbstractDomain}

An outer constructor for the `BucketModel`, which enforces the
constraints:
1. The bucket model domain is of type <: ClimaLand.Domains.AbstractDomain
2. Using an albedo read from a lat/lon file requires a global run.
"""
function BucketModel(;
    parameters::BucketModelParameters{FT, PSE},
    domain::D,
    atmosphere::ATM,
    radiation::RAD,
) where {FT, PSE, ATM, RAD, D <: ClimaLand.Domains.AbstractDomain}
    if parameters.albedo isa Union{BulkAlbedoStatic, BulkAlbedoTemporal}
        typeof(domain) <: SphericalShell ? nothing :
        error("Using an albedo map requires a global run.")
    end

    args = (parameters, atmosphere, radiation, domain)
    BucketModel{FT, typeof.(args)...}(args...)
end



prognostic_types(::BucketModel{FT}) where {FT} = (FT, FT, FT, FT)
prognostic_vars(::BucketModel) = (:W, :T, :Ws, :σS)
prognostic_domain_names(::BucketModel) =
    (:surface, :subsurface, :surface, :surface)

auxiliary_types(::BucketModel{FT}) where {FT} = (
    FT,
    NamedTuple{(:lhf, :shf, :vapor_flux, :r_ae), Tuple{FT, FT, FT, FT}},
    FT,
    FT,
    FT,
    FT,
)
auxiliary_vars(::BucketModel) =
    (:q_sfc, :turbulent_fluxes, :R_n, :T_sfc, :α_sfc, :ρ_sfc)
auxiliary_domain_names(::BucketModel) =
    (:surface, :surface, :surface, :surface, :surface, :surface, :surface)

"""
    ClimaLand.make_set_initial_cache(model::BucketModel{FT}) where{FT}

Returns the set_initial_cache! function, which updates the auxiliary
state `p` in place with the initial values corresponding to Y(t=t0) = Y0.

In this case, we also use this method to update the initial values for the
spatially varying parameter fields, read in from data files.
"""
function ClimaLand.make_set_initial_cache(model::BucketModel)
    update_cache! = make_update_cache(model)
    function set_initial_cache!(p, Y0, t0)
        update_cache!(p, Y0, t0)
    end
    return set_initial_cache!
end

"""
    make_compute_exp_tendency(model::BucketModel{FT}) where {FT}

Creates the compute_exp_tendency! function for the bucket model.
"""
function make_compute_exp_tendency(model::BucketModel{FT}) where {FT}
    function compute_exp_tendency!(dY, Y, p, t)
        (; κ_soil, ρc_soil, σS_c, W_f) = model.parameters

        #Currently, the entire surface is assumed to be
        # snow covered entirely or not at all.
        snow_cover_fraction = heaviside.(Y.bucket.σS)

        # In this case, there is just one set of surface fluxes to compute.
        # Since q_sfc itself has already been modified to account for
        # snow covered regions, and since the surface temperature is
        # assumed to be the same for snow and underlying land,
        # the returned fluxes are computed correctly for the cell
        # regardless of snow-coverage.

        # The below must be adjusted to compute F_sfc over snow and over soil
        # if we want the snow cover fraction to be intermediate between 0 and 1.
        (; turbulent_fluxes, R_n) = p.bucket
        F_sfc = @. (turbulent_fluxes.shf .+ turbulent_fluxes.lhf + R_n) # Eqn (21)
        _T_freeze = LP.T_freeze(model.parameters.earth_param_set)
        _LH_f0 = LP.LH_f0(model.parameters.earth_param_set)
        _ρ_liq = LP.ρ_cloud_liq(model.parameters.earth_param_set)
        _ρLH_f0 = _ρ_liq * _LH_f0 # Latent heat per unit volume.
        # partition energy fluxes for snow covered area
        partitioned_fluxes =
            partition_snow_surface_fluxes.(
                Y.bucket.σS,
                p.bucket.T_sfc,
                model.parameters.τc,
                snow_cover_fraction,
                turbulent_fluxes.vapor_flux,
                F_sfc,
                _ρLH_f0,
                _T_freeze,
            )
        (; F_melt, F_into_snow, G_under_snow) = partitioned_fluxes
        G = @. F_sfc * (1 - snow_cover_fraction) +
           G_under_snow * snow_cover_fraction # Equation 22
        # Temperature profile of soil.
        gradc2f = ClimaCore.Operators.GradientC2F()
        divf2c = ClimaCore.Operators.DivergenceF2C(
            top = ClimaCore.Operators.SetValue(ClimaCore.Geometry.WVector.(G)),
            bottom = ClimaCore.Operators.SetValue(
                ClimaCore.Geometry.WVector.(FT(0.0)),
            ),
        )

        @. dY.bucket.T = -1 / ρc_soil * (divf2c(-κ_soil * gradc2f(Y.bucket.T))) # Simple heat equation, Eq 10

        # Partition water fluxes
        liquid_precip = liquid_precipitation(model.atmos, p, t) # always negative
        snow_precip = snow_precipitation(model.atmos, p, t) # always negative
        # F_melt is negative as it is a downward welling flux warming the snow
        snow_melt = @. F_melt / _ρLH_f0 # defined after Equation (22)

        infiltration = @. infiltration_at_point(
            Y.bucket.W,
            snow_cover_fraction * snow_melt,
            liquid_precip,
            (1 - snow_cover_fraction) * turbulent_fluxes.vapor_flux,
            W_f,
        ) # Equation (2) of the text.

        # Positive infiltration -> net (negative) flux into soil
        @. dY.bucket.W = -infiltration # Equation (2) of the text.

        dY.bucket.Ws = @. -(
            liquid_precip +
            snow_cover_fraction * snow_melt +
            (1 - snow_cover_fraction) * turbulent_fluxes.vapor_flux -
            infiltration
        ) # Equation (3) of the text

        dY.bucket.σS = @. -(
            snow_precip + snow_cover_fraction * turbulent_fluxes.vapor_flux -
            snow_cover_fraction * snow_melt
        ) # Equation (6)
    end
    return compute_exp_tendency!
end

"""
    make_update_aux(model::BucketModel{FT}) where {FT}

Creates the update_aux! function for the BucketModel.
"""
function make_update_aux(model::BucketModel{FT}) where {FT}
    function update_aux!(p, Y, t)
        p.bucket.T_sfc .= ClimaLand.Domains.top_center_to_surface(Y.bucket.T)
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
                Ref(
                    LP.thermodynamic_parameters(
                        model.parameters.earth_param_set,
                    ),
                ),
            )

        # Compute turbulent surface fluxes
        p.bucket.turbulent_fluxes .=
            turbulent_fluxes(model.atmos, model, Y, p, t)

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
of snow and of the bareground surface, based on the snow water equivalent `S` relative to
the parameter `S_c`. The linear interpolation is taken from Lague et al 2019.

Note that if our snow_cover_fraction function was smoothly varying, the albedo
would simply be σα_snow + (1-σ)α_bareground. Since we cannot support snow cover
fractions that are not a heaviside function, we have a small inconsistency
for 0 < σS < eps(FT) where the snow cover fraction is zero, but there is a small
contribution of snow to the albedo.
"""
function next_albedo(
    model_albedo::Union{BulkAlbedoFunction{FT}, BulkAlbedoStatic{FT}},
    parameters,
    Y,
    p,
    t,
) where {FT}
    (; α_snow, α_bareground) = model_albedo
    (; σS_c) = parameters
    σS = max.(Y.bucket.σS, FT(0))# if σS is small and negative, clip to zero for albedo estimation
    return @. (
        (1 - σS / (σS + σS_c)) * α_bareground + σS / (σS + σS_c) * α_snow
    )
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
    # Read next data fields if initializing or next date is closest to current time
    # This maintains `all_dates[date_idx]` <= `sim_date` < `all_dates[date_idx + 1]`
    if t == sim_info.t_start ||
       sim_date >= to_datetime(next_date_in_file(model_albedo.albedo_info))
        read_data_fields!(model_albedo.albedo_info, sim_date, axes(Y.bucket.W))
    end
    # Interpolate data value to current time
    return get_data_at_date(
        model_albedo.albedo_info,
        axes(Y.bucket.W),
        "sw_alb",
        sim_date,
    )
end

function ClimaLand.get_drivers(model::BucketModel)
    return (model.atmos, model.radiation)
end

include("./bucket_parameterizations.jl")

end

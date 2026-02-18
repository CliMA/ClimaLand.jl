import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import Interpolations: Constant
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.FileReaders: NCFileReader, read
using DocStringExtensions

export prescribed_lai_era5,
    prescribed_lai_modis,
    PrescribedAreaIndices,
    update_biomass!,
    ZhouOptimalLAIModel



"""
     prescribed_lai_era5(era5_lai_ncdata_path,
                         era5_lai_cover_ncdata_path,
                         surface_space,
                         start_date,
                         earth_param_set;
                         time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
                         regridder_type = :InterpolationsRegridder,
                         interpolation_method = Interpolations.Constant(),)

A helper function which constructs the TimeVaryingInput object for Leaf Area Index, from a
file path pointing to the ERA5 LAI data in a netcdf file, a file path pointing to the ERA5
LAI cover data in a netcdf file, the surface_space, the start date, and the earth_param_set.

This currently one works when a single file is passed for both the era5 lai and era5 lai cover data.

The ClimaLand default is to use nearest neighbor interpolation, but
linear interpolation is supported
by passing interpolation_method = Interpolations.Linear().
"""
function prescribed_lai_era5(
    era5_lai_ncdata_path,
    era5_lai_cover_ncdata_path,
    surface_space,
    start_date;
    time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
    regridder_type = :InterpolationsRegridder,
    interpolation_method = Interpolations.Constant(),
)
    hvc_ds = NCFileReader(era5_lai_cover_ncdata_path, "cvh")
    lvc_ds = NCFileReader(era5_lai_cover_ncdata_path, "cvl")
    hv_cover = read(hvc_ds)
    lv_cover = read(lvc_ds)
    close(hvc_ds)
    close(lvc_ds)
    compose_function = let hv_cover = hv_cover, lv_cover = lv_cover
        (lai_hv, lai_lv) -> lai_hv .* hv_cover .+ lai_lv .* lv_cover
    end
    return TimeVaryingInput(
        era5_lai_ncdata_path,
        ["lai_hv", "lai_lv"],
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        method = time_interpolation_method,
        compose_function = compose_function,
    )
end

"""
     prescribed_lai_modis(surface_space,
                          start_date
                          stop_date,
                          earth_param_set;
                          time_interpolation_method = LinearInterpolation(),
                          regridder_type = :InterpolationsRegridder,
                          interpolation_method = Interpolations.Constant(),
                          modis_lai_ncdata_path = nothing,
                          context = ClimaComms.context(surface_space))

A helper function which constructs the TimeVaryingInput object for Leaf Area
Index using MODIS LAI data; requires the
surface_space, the start and stop dates, and the earth_param_set.

The ClimaLand default is to use nearest neighbor interpolation, but
linear interpolation is supported
by passing interpolation_method = Interpolations.Linear().

If `modis_lai_ncdata_path` is provided, it will be used directly.
Otherwise, the path will be inferred from the start and stop dates.
"""
function prescribed_lai_modis(
    surface_space,
    start_date,
    stop_date;
    time_interpolation_method = LinearInterpolation(),
    regridder_type = :InterpolationsRegridder,
    interpolation_method = Interpolations.Constant(),
    modis_lai_ncdata_path = nothing,
    context = ClimaComms.context(surface_space),
)
    modis_lai_ncdata_path =
        isnothing(modis_lai_ncdata_path) ?
        ClimaLand.Artifacts.modis_lai_multiyear_paths(;
            context,
            start_date,
            stop_date,
        ) : modis_lai_ncdata_path
    return TimeVaryingInput(
        modis_lai_ncdata_path,
        ["lai"],
        surface_space;
        start_date,
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        method = time_interpolation_method,
    )
end


"""
    AbstractBiomassModel{FT} <: AbstractCanopyComponent{FT}

An abstract type for modeling the biomass (above ground - LAI, SAI, canopy
height) and below ground (rooting depth, RAI).
"""
abstract type AbstractBiomassModel{FT} <: AbstractCanopyComponent{FT} end

ClimaLand.name(::AbstractBiomassModel) = :biomass

abstract type AbstractAreaIndexModel{FT <: AbstractFloat} end

"""
   PrescribedAreaIndices{FT <:AbstractFloat, F <: AbstractTimeVaryingInput}

A struct containing the area indices of the plants at a specific site;
LAI varies in time, while SAI and RAI are fixed in space and time.

$(DocStringExtensions.FIELDS)
"""
struct PrescribedAreaIndices{
    FT <: AbstractFloat,
    F <: AbstractTimeVaryingInput,
} <: AbstractAreaIndexModel{FT}
    "A function of simulation time `t` giving the leaf area index (LAI; m2/m2)"
    LAI::F
    "The constant stem area index (SAI; m2/m2)"
    SAI::FT
    "The constant root area index (RAI; m2/m2)"
    RAI::FT
end

"""
    PrescribedAreaIndices{FT}(
        LAI::AbstractTimeVaryingInput,
        SAI::FT,
        RAI::FT,
    ) where {FT <: AbstractFloat}

An outer constructor for setting the PrescribedAreaIndices given
LAI, SAI, and RAI.
"""
function PrescribedAreaIndices{FT}(
    LAI::AbstractTimeVaryingInput,
    SAI::FT,
    RAI::FT,
) where {FT <: AbstractFloat}
    PrescribedAreaIndices{FT, typeof(LAI)}(LAI, SAI, RAI)
end

"""
    lai_consistency_check(
        n_stem::Int64,
        n_leaf::Int64,
        area_index::NamedTuple{(:root, :stem, :leaf), Tuple{FT, FT, FT}},
    ) where {FT}

Carries out consistency checks using the area indices supplied and the number of
stem elements being modeled.

Note that it is possible to have a plant with no stem compartments
but with leaf compartments, and that a plant must have leaf compartments
(even if LAI = 0).

Specifically, this checks that:
1. n_leaf > 0
2. if LAI is nonzero or SAI is nonzero, RAI must be nonzero.
3. if SAI > 0, n_stem must be > 0 for consistency. If SAI == 0, n_stem must
be zero.
"""
function lai_consistency_check(
    n_stem::Int64,
    n_leaf::Int64,
    area_index::NamedTuple{(:root, :stem, :leaf), Tuple{FT, FT, FT}},
) where {FT}
    @assert n_leaf > 0
    if area_index.leaf > eps(FT) || area_index.stem > eps(FT)
        @assert area_index.root > eps(FT)
    end
    # If there SAI is zero, there should be no stem compartment
    if area_index.stem < eps(FT)
        @assert n_stem == FT(0)
    else
        # if SAI is > 0, n_stem should be > 0 for consistency
        @assert n_stem > 0
    end

end

"""
    struct PrescribedBiomassModel{FT, PSAI <: PrescribedAreaIndices, RDTH <: Union{FT, ClimaCore.Fields.Field}} <: AbstractBiomassModel{FT}

A prescribed biomass model where LAI, SAI, RAI, rooting depth, and height are prescribed.

In  global run with patches
of bare soil, you can "turn off" the canopy model (to get zero root extraction, zero absorption and
emission, zero transpiration and sensible heat flux from the canopy), by setting:
- n_leaf = 1
- n_stem = 0
- LAI = SAI = RAI = 0.

If run with PlantHydraulics, this must be consistent with the plant hydraulics model:a
 plant model can have leaves but no stem, but not vice versa. If n_stem = 0, SAI must be zero.

$(DocStringExtensions.FIELDS)
"""
struct PrescribedBiomassModel{
    FT,
    PSAI <: PrescribedAreaIndices{FT},
    RDTH <: Union{FT, ClimaCore.Fields.Field},
    HTH <: Union{FT, ClimaCore.Fields.Field},
} <: AbstractBiomassModel{FT}
    "The plant area index model for LAI, SAI, RAI"
    plant_area_index::PSAI
    "Rooting depth parameter (m) - a characteristic depth below which 1/e of the root mass lies"
    rooting_depth::RDTH
    "Canopy height (m) - can be scalar (uniform) or spatially-varying Field"
    height::HTH
end

"""
    PrescribedBiomassModel{FT}(;LAI::AbstractTimeVaryingInput,
                                SAI::FT,
                                RAI::FT,
                                rooting_depth,
                                height) where {FT}

An outer constructor to help set up the PrescribedBiomassModel from 
LAI, SAI, and RAI directly, instead of requiring the user to make the
area index object first; rooting_depth and height are also required.

Height can be either:
- A scalar FT value (uniform height across domain)
- A ClimaCore.Fields.Field (spatially-varying height)
"""
function PrescribedBiomassModel{FT}(;
    LAI::AbstractTimeVaryingInput,
    SAI::FT,
    RAI::FT,
    rooting_depth,
    height,
) where {FT}
    plant_area_index = PrescribedAreaIndices{FT}(LAI, SAI, RAI)
    args = (plant_area_index, rooting_depth, height)
    PrescribedBiomassModel{FT, typeof.(args)...}(args...)
end

ClimaLand.auxiliary_vars(model::PrescribedBiomassModel) = (:area_index,)
ClimaLand.auxiliary_types(model::PrescribedBiomassModel{FT}) where {FT} =
    (NamedTuple{(:root, :stem, :leaf), Tuple{FT, FT, FT}},)
ClimaLand.auxiliary_domain_names(::PrescribedBiomassModel) = (:surface,)

function clip(x::FT, threshold::FT) where {FT}
    x > threshold ? x : FT(0)
end

"""
    update_biomass!(
        p,
        Y,
        t,
        component::PrescribedBiomassModel{FT},
        canopy,
    ) where {FT}

Sets the area indices pertaining to their values at time t.

Note that we clip all values of LAI below 0.05 to zero.
This is because we currently run into issues when LAI is
of order eps(FT) in the SW radiation code.
Please see Issue #644
or PR #645 for details.
For now, this clipping is similar to what CLM and NOAH MP do.
"""

function update_biomass!(
    p,
    Y,
    t,
    component::PrescribedBiomassModel{FT},
    canopy,
) where {FT}
    (; LAI, SAI, RAI) = component.plant_area_index
    evaluate!(p.canopy.biomass.area_index.leaf, LAI, t)
    p.canopy.biomass.area_index.leaf .=
        clip.(p.canopy.biomass.area_index.leaf, FT(0.05))
    @. p.canopy.biomass.area_index.stem = SAI
    @. p.canopy.biomass.area_index.root = RAI
end

"""
    root_distribution(z::FT, rooting_depth::FT)

Computes value of rooting probability density function at `z`.

The rooting probability density function is derived from the
cumulative distribution function F(z) = 1 - β^(100z), which is described
by Equation 2.23 of
Bonan, "Climate Change and Terrestrial Ecosystem Modeling", 2019 Cambridge University Press.
This probability distribution function is equivalent to the derivative of the
cumulative distribution function with respect to z,
where `rooting_depth` replaces (-1)/(100ln(β)) and z is expected to be negative.
"""
function root_distribution(z::FT, rooting_depth::FT) where {FT <: AbstractFloat}
    return (1 / rooting_depth) * exp(z / rooting_depth) # 1/m
end

#####################################################################
# ZhouOptimalLAIModel - Optimal LAI model based on Zhou et al. (2025)
#####################################################################

"""
    ZhouOptimalLAIModel{FT, OLPT <: OptimalLAIParameters{FT}, GD, RDTH, HTH} <: AbstractBiomassModel{FT}

An implementation of the optimal LAI model from Zhou et al. (2025) as a biomass model.

This model computes LAI dynamically based on optimality principles, balancing energy and
water constraints. LAI is stored in `p.canopy.biomass.area_index.leaf`, consistent with
`PrescribedBiomassModel`.

# Fields
- `parameters`: Required parameters for the optimal LAI model
- `optimal_lai_inputs`: NamedTuple with spatially varying GSL (growing season length in days),
  A0_annual (annual potential GPP in mol CO2 m^-2 yr^-1), precip_annual (mol H2O m^-2 yr^-1),
  vpd_gs (Pa), lai_init (initial LAI from MODIS), and f0 (fraction of precip for transpiration),
  typically created using `optimal_lai_initial_conditions`.
- `SAI`: Prescribed stem area index (m^2 m^-2)
- `RAI`: Prescribed root area index (m^2 m^-2)
- `rooting_depth`: Rooting depth parameter (m) - a characteristic depth below which 1/e of the root mass lies
- `height`: Canopy height (m) - can be scalar (uniform) or spatially-varying Field

# References
Zhou et al. (2025) "A General Model for the Seasonal to Decadal Dynamics of Leaf Area"
Global Change Biology. https://onlinelibrary.wiley.com/doi/pdf/10.1111/gcb.70125
"""
struct ZhouOptimalLAIModel{
    FT,
    OLPT <: OptimalLAIParameters{FT},
    GD,
    RDTH <: Union{FT, ClimaCore.Fields.Field},
    HTH <: Union{FT, ClimaCore.Fields.Field},
} <: AbstractBiomassModel{FT}
    "Required parameters for the optimal LAI model"
    parameters::OLPT
    "Spatially varying initial conditions (GSL, A0_annual, precip_annual, vpd_gs, lai_init, f0)"
    optimal_lai_inputs::GD
    "Prescribed stem area index (m^2 m^-2)"
    SAI::FT
    "Prescribed root area index (m^2 m^-2)"
    RAI::FT
    "Rooting depth parameter (m)"
    rooting_depth::RDTH
    "Canopy height (m) - can be scalar (uniform) or spatially-varying Field"
    height::HTH
end

Base.eltype(::ZhouOptimalLAIModel{FT}) where {FT} = FT

"""
    ZhouOptimalLAIModel{FT}(
        parameters::OptimalLAIParameters{FT},
        optimal_lai_inputs;
        SAI::FT,
        RAI::FT,
        rooting_depth,
        height,
    ) where {FT <: AbstractFloat}

Outer constructor for the ZhouOptimalLAIModel struct.

# Arguments
- `parameters`: OptimalLAIParameters for the model
- `optimal_lai_inputs`: NamedTuple with spatially varying GSL, A0_annual, precip_annual, vpd_gs, lai_init, f0 fields,
  typically created using `optimal_lai_initial_conditions`.
- `SAI`: Prescribed stem area index (m^2 m^-2)
- `RAI`: Prescribed root area index (m^2 m^-2)
- `rooting_depth`: Rooting depth parameter (m)
- `height`: Canopy height (m) - can be scalar or spatially-varying Field
"""
function ZhouOptimalLAIModel{FT}(
    parameters::OptimalLAIParameters{FT},
    optimal_lai_inputs;
    SAI::FT,
    RAI::FT,
    rooting_depth,
    height,
) where {FT <: AbstractFloat}
    return ZhouOptimalLAIModel{
        FT,
        typeof(parameters),
        typeof(optimal_lai_inputs),
        typeof(rooting_depth),
        typeof(height),
    }(
        parameters,
        optimal_lai_inputs,
        SAI,
        RAI,
        rooting_depth,
        height,
    )
end

"""
    ClimaLand.auxiliary_vars(model::ZhouOptimalLAIModel)
    ClimaLand.auxiliary_types(model::ZhouOptimalLAIModel)
    ClimaLand.auxiliary_domain_names(model::ZhouOptimalLAIModel)

Defines the auxiliary variables for the ZhouOptimalLAIModel:
- `area_index`: NamedTuple{(:root, :stem, :leaf)} containing area indices (m^2 m^-2)
- `A0_daily`: daily potential GPP from previous day (mol CO2 m^-2 day^-1), with moisture stress factor beta
- `A0_annual`: annual potential GPP from previous year (mol CO2 m^-2 yr^-1), with actual moisture stress factor beta
- `A0_daily_acc`: accumulator for current day's potential GPP (mol CO2 m^-2 day^-1)
- `A0_annual_acc`: accumulator for current year's potential GPP (mol CO2 m^-2 yr^-1), sum of daily values at noon
- `GSL`: growing season length (days), spatially varying
- `precip_annual`: mean annual precipitation (mol H2O m^-2 yr^-1), for water limitation in LAI_max
- `vpd_gs`: mean VPD during growing season (Pa), for water limitation WUE factor in LAI_max
- `f0`: spatially varying fraction of precipitation for transpiration (dimensionless), from Zhou et al.
"""
ClimaLand.auxiliary_vars(model::ZhouOptimalLAIModel) = (
    :area_index,
    :A0_daily,
    :A0_annual,
    :A0_daily_acc,
    :A0_annual_acc,
    :days_since_reset,
    :GSL,
    :precip_annual,
    :vpd_gs,
    :f0,
)
ClimaLand.auxiliary_types(model::ZhouOptimalLAIModel{FT}) where {FT} = (
    NamedTuple{(:root, :stem, :leaf), Tuple{FT, FT, FT}},
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
)
ClimaLand.auxiliary_domain_names(::ZhouOptimalLAIModel) = (
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
    :surface,
)

"""
    update_biomass!(
        p,
        Y,
        t,
        component::ZhouOptimalLAIModel{FT},
        canopy,
    ) where {FT}

Sets the SAI and RAI to their prescribed constant values.

Note: LAI is updated via the optimal LAI callback (make_OptimalLAI_callback) at local noon,
not in this function. This function only handles the constant area indices (SAI and RAI).
"""
function update_biomass!(
    p,
    Y,
    t,
    component::ZhouOptimalLAIModel{FT},
    canopy,
) where {FT}
    (; SAI, RAI) = component
    @. p.canopy.biomass.area_index.stem = SAI
    @. p.canopy.biomass.area_index.root = RAI
    # LAI is updated via the callback, not here
    # Apply clipping to LAI (same as PrescribedBiomassModel)
    p.canopy.biomass.area_index.leaf .=
        clip.(p.canopy.biomass.area_index.leaf, FT(0.05))
end

"""
    set_historical_cache!(p, Y0, model::ZhouOptimalLAIModel, canopy; A0_daily)

The optimal LAI model requires initialization of LAI and A0 values before the simulation.

GSL, A0_annual, precip_annual, vpd_gs, lai_init, and f0 are taken from `model.optimal_lai_inputs`, which
contains spatially varying fields typically created using `optimal_lai_initial_conditions`.

LAI is initialized from MODIS satellite observations (`lai_init`) rather than computing
equilibrium from model equations. This provides a realistic starting point that reduces
spin-up time and matches observed vegetation patterns.

# Arguments
- `A0_daily`: Initial daily potential GPP (mol CO2 m^-2 day^-1), default 0.5.

# Notes
- precip_annual is in mol H2O m^-2 yr^-1 (converted from ERA5 kg m^-2 s^-1 in optimal_lai_inputs.jl)
  to match Zhou et al. (2025) formulation for water limitation.
- vpd_gs is the mean VPD during growing season (Pa), used in the WUE factor for water limitation.
- lai_init comes from MODIS first timestep, providing spatially-varying initial conditions.
- f0 is the spatially varying fraction of precipitation for transpiration from Zhou et al.
"""
function set_historical_cache!(
    p,
    Y0,
    model::ZhouOptimalLAIModel,
    canopy;
    A0_daily = zero(eltype(model.parameters)),
)
    parameters = model.parameters
    optimal_lai_inputs = model.optimal_lai_inputs
    FT = eltype(parameters)

    # Get spatially varying data from optimal_lai_inputs
    GSL = optimal_lai_inputs.GSL
    A0_annual = optimal_lai_inputs.A0_annual
    precip_annual = optimal_lai_inputs.precip_annual  # in mol H2O m^-2 yr^-1
    vpd_gs = optimal_lai_inputs.vpd_gs  # in Pa
    lai_init = optimal_lai_inputs.lai_init  # MODIS first timestep
    f0 = optimal_lai_inputs.f0  # spatially varying f0 from Zhou et al.

    L = p.canopy.biomass.area_index.leaf

    # Initialize LAI from MODIS satellite observations
    # This provides realistic spatially-varying initial conditions that reduce spin-up time
    L .= lai_init

    # Initialize A0 variables (supports both scalar and Field inputs via .=)
    p.canopy.biomass.A0_daily .= A0_daily
    p.canopy.biomass.A0_annual .= A0_annual
    p.canopy.biomass.A0_daily_acc .= FT(0)
    p.canopy.biomass.A0_annual_acc .= FT(0)
    p.canopy.biomass.days_since_reset .= FT(0)

    # Store GSL in the cache (spatially varying field)
    p.canopy.biomass.GSL .= GSL

    # Store precip_annual (already in mol H2O m^-2 yr^-1 from optimal_lai_inputs.nc)
    p.canopy.biomass.precip_annual .= precip_annual

    # Store vpd_gs (already in Pa)
    p.canopy.biomass.vpd_gs .= vpd_gs

    # Store f0 (spatially varying fraction of precipitation for transpiration)
    p.canopy.biomass.f0 .= f0
end

"""
    get_model_callbacks(component::ZhouOptimalLAIModel, canopy; t0, Δt)

Creates the optimal LAI callback and returns it as a single element tuple of model callbacks.

# Notes
- Daily A0 is computed with moisture stress factor beta - drives L_steady
- Annual A0 is computed with actual moisture stress factor beta - used for LAI_max
- Water limitation enters LAI_max through the f0*P/A0 * (ca(1-chi))/(1.6*D) term (Equation 11)
- GSL (Growing Season Length) is read from p.canopy.biomass.GSL, which supports spatially
  varying values initialized via set_historical_cache!.
"""
function get_model_callbacks(
    component::ZhouOptimalLAIModel{FT},
    canopy;
    t0,
    Δt,
) where {FT}
    lai_cb = make_OptimalLAI_callback(FT, t0, Δt, canopy)
    return (lai_cb,)
end

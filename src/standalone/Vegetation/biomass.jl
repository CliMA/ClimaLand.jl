import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, LinearInterpolation, PeriodicCalendar
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import Interpolations: Constant
import ClimaUtilities.Regridders: InterpolationsRegridder
import ClimaUtilities.FileReaders: NCFileReader, read
using DocStringExtensions

export prescribed_lai_era5,
    prescribed_lai_modis,
    prescribed_climatological_lai_modis,
    PrescribedAreaIndices,
    update_biomass!,
    ZhouOptimalLAIModel,
    mask_biomass!

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
by passing `interpolation_method = Interpolations.Linear()`.
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
by passing `interpolation_method = Interpolations.Linear()`.

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
     prescribed_climatological_lai_modis(surface_space,
                                         time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
                                         regridder_type = :InterpolationsRegridder,
                                         interpolation_method = Interpolations.Constant(),
                                         context = ClimaComms.context(surface_space))

A helper function which constructs the TimeVaryingInput object for Leaf Area
Index using MODIS climatological LAI data; requires the
surface_space.

The ClimaLand default is to use nearest neighbor interpolation, but
linear interpolation is supported
by passing interpolation_method = Interpolations.Linear().
"""
function prescribed_climatological_lai_modis(
    surface_space,
    time_interpolation_method = LinearInterpolation(PeriodicCalendar()),
    regridder_type = :InterpolationsRegridder,
    interpolation_method = Interpolations.Constant(),
    context = ClimaComms.context(surface_space),
)
    modis_lai_ncdata_path =
        ClimaLand.Artifacts.modis_lai_climatology_data_path(; context)
    return TimeVaryingInput(
        modis_lai_ncdata_path,
        ["lai"],
        surface_space;
        regridder_type,
        regridder_kwargs = (; interpolation_method),
        method = time_interpolation_method,
    )
end


"""
     modis_max_lai(surface_space,
                   regridder_type = :InterpolationsRegridder,
                   interpolation_method = Interpolations.Constant(),
                   context = ClimaComms.context(surface_space))

A helper function which constructs the SpaceVaryingInput object for the maximum
Leaf Area Index using MODIS LAI data; requires the surface_space.

The ClimaLand default is to use nearest neighbor interpolation, but
linear interpolation is supported by passing
`interpolation_method = Interpolations.Linear()`.
"""
function modis_max_lai(
    surface_space,
    regridder_type = :InterpolationsRegridder,
    interpolation_method = Interpolations.Constant(),
    context = ClimaComms.context(surface_space),
)
    modis_max_lai_ncdata_path =
        ClimaLand.Artifacts.modis_max_lai_data_path(; context)
    return SpaceVaryingInput(
        modis_max_lai_ncdata_path,
        ["lai"],
        surface_space;
        regridder_type,
        regridder_kwargs = (; interpolation_method),
    )
end

"""
    AbstractBiomassModel{FT} <: AbstractCanopyComponent{FT}

An abstract type for modeling the biomass (above ground - LAI, SAI, canopy
height) and below ground (rooting depth, RAI).
"""
abstract type AbstractBiomassModel{FT} <: AbstractCanopyComponent{FT} end

ClimaLand.name(::AbstractBiomassModel) = :biomass

abstract type AbstractAreaIndexModel end

"""
   PrescribedAreaIndices{FS <: Union{AbstractFloat, ClimaCore.Fields.Field}, F <: AbstractTimeVaryingInput}

A struct containing the area indices of the plants at a specific site;
LAI varies in time, while SAI and RAI are fixed in time and can either be a
scalar (spatially uniform) or a ClimaCore Field (spatially varying).

$(DocStringExtensions.FIELDS)
"""
struct PrescribedAreaIndices{
    FS <: Union{AbstractFloat, ClimaCore.Fields.Field},
    F <: AbstractTimeVaryingInput,
} <: AbstractAreaIndexModel
    "A function of simulation time `t` giving the leaf area index (LAI; m2/m2)"
    LAI::F
    "The constant-in-time stem area index (SAI; m2/m2), scalar or Field"
    SAI::FS
    "The constant-in-time root area index (RAI; m2/m2), scalar or Field"
    RAI::FS
end

"""
    PrescribedAreaIndices(
        LAI::AbstractTimeVaryingInput,
        SAI,
        RAI,
    )

An outer constructor for setting the PrescribedAreaIndices given LAI, SAI, and
RAI. SAI and RAI may be scalars or ClimaCore Fields.
"""
function PrescribedAreaIndices(LAI::AbstractTimeVaryingInput, SAI, RAI)
    PrescribedAreaIndices{typeof(SAI), typeof(LAI)}(LAI, SAI, RAI)
end

"""
    struct PrescribedBiomassModel{FT, PSAI <: PrescribedAreaIndices, RDTH <: Union{FT, ClimaCore.Fields.Field}} <: AbstractBiomassModel{FT}

A prescribed biomass model where LAI, SAI, RAI, rooting depth, and height are prescribed.

In  global run with patches
of bare soil, you can "turn off" the canopy model (to get zero root extraction, zero absorption and
emission, zero transpiration and sensible heat flux from the canopy), by setting:
- LAI = SAI = RAI = 0.
$(DocStringExtensions.FIELDS)
"""
struct PrescribedBiomassModel{
    FT,
    PSAI <: PrescribedAreaIndices,
    RDTH <: Union{FT, ClimaCore.Fields.Field},
    HTH <: Union{FT, ClimaCore.Fields.Field},
} <: AbstractBiomassModel{FT}
    "The plant area index model for LAI, SAI, RAI"
    plant_area_index::PSAI
    "Rooting depth parameter (m) - a characteristic depth below which 1/e of the root mass lies"
    rooting_depth::RDTH
    "Canopy height (m) - can be scalar (uniform) or spatially-varying Field"
    height::HTH
    function PrescribedBiomassModel{FT, PSAI, RDTH, HTH}(
        plant_area_index,
        rooting_depth,
        height,
    ) where {FT, PSAI, RDTH, HTH}
        new{FT, PSAI, RDTH, HTH}(plant_area_index, rooting_depth, height)
    end
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
    SAI,
    RAI,
    rooting_depth,
    height,
) where {FT}
    plant_area_index = PrescribedAreaIndices(LAI, SAI, RAI)
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
    mask_biomass!(p, Val(canopy.boundary_conditions.prognostic_land_components))
end

"""
    mask_biomass!(p, prognostic_land_components)

Default method of setting LAI/RAI/SAI to zero where there
cannot be canopy; does nothing.

Currently, this is only does something when a lake model
is included in integrated models, as we cannot have 
vegetation over a lake, and the lake masks may not be consistent with
the biomass model.
"""
mask_biomass!(p, prognostic_land_components) = nothing

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
water constraints. LAI is prognostic, in `Y.canopy.biomass.LAI`, and is mirrored into
`p.canopy.biomass.area_index.leaf` for the rest of the canopy, consistent with
`PrescribedBiomassModel`. Its initial state, and that of the trailing totals it depends
on, is set by `ClimaLand.Simulations.set_canopy_biomass_initial_conditions!`.

# Fields
- `parameters`: Required parameters for the optimal LAI model
- `optimal_lai_inputs`: NamedTuple with the spatially varying inputs of the LAI_max and
  steady-state-LAI formulas: GSL (growing season length in days), vpd_gs (growing-season
  mean VPD in Pa), and f0 (fraction of precip for transpiration). Typically created using
  `optimal_lai_initial_conditions`, which reads them alongside the initial conditions.
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
    FS <: Union{FT, ClimaCore.Fields.Field},
    RDTH <: Union{FT, ClimaCore.Fields.Field},
    HTH <: Union{FT, ClimaCore.Fields.Field},
} <: AbstractBiomassModel{FT}
    "Required parameters for the optimal LAI model"
    parameters::OLPT
    "Spatially varying inputs of the LAI formulas (GSL, vpd_gs, f0)"
    optimal_lai_inputs::GD
    "Prescribed stem area index (m^2 m^-2), scalar or Field"
    SAI::FS
    "Prescribed root area index (m^2 m^-2), scalar or Field"
    RAI::FS
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
        SAI,
        RAI,
        rooting_depth,
        height,
    ) where {FT <: AbstractFloat}

Outer constructor for the ZhouOptimalLAIModel struct.

# Arguments
- `parameters`: OptimalLAIParameters for the model
- `optimal_lai_inputs`: NamedTuple with the spatially varying GSL, vpd_gs and f0 fields,
  typically created using `optimal_lai_initial_conditions`.
- `SAI`: Prescribed stem area index (m^2 m^-2); scalar or spatially-varying Field
- `RAI`: Prescribed root area index (m^2 m^-2); scalar or spatially-varying Field
- `rooting_depth`: Rooting depth parameter (m)
- `height`: Canopy height (m) - can be scalar or spatially-varying Field
"""
function ZhouOptimalLAIModel{FT}(
    parameters::OptimalLAIParameters{FT},
    optimal_lai_inputs;
    SAI,
    RAI,
    rooting_depth,
    height,
) where {FT <: AbstractFloat}
    return ZhouOptimalLAIModel{
        FT,
        typeof(parameters),
        typeof(optimal_lai_inputs),
        typeof(SAI),
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
- `A0_inst`: instantaneous potential GPP (mol CO2 m^-2 s^-1), computed once per tendency evaluation and shared by the daily and annual A0 running sums
- `GSL`: growing season length (days), spatially varying
- `vpd_gs`: mean VPD during growing season (Pa), for water limitation WUE factor in LAI_max
- `f0`: spatially varying fraction of precipitation for transpiration (dimensionless), from Zhou et al.
"""
ClimaLand.auxiliary_vars(model::ZhouOptimalLAIModel) =
    (:area_index, :A0_inst, :GSL, :vpd_gs, :f0, :f0_base)
ClimaLand.auxiliary_types(model::ZhouOptimalLAIModel{FT}) where {FT} =
    (NamedTuple{(:root, :stem, :leaf), Tuple{FT, FT, FT}}, FT, FT, FT, FT, FT)
ClimaLand.auxiliary_domain_names(::ZhouOptimalLAIModel) =
    (:surface, :surface, :surface, :surface, :surface, :surface)

# The optimal-LAI model's prognostic state is nine time-integrated variables in `Y`,
# all advanced smoothly by the time-stepper (no callback, checkpoint/restart-safe):
#   `A0_daily`      1-day potential-GPP total (RunningSum, τ = 1 day),
#   `A0_annual`     smoothed 1-year potential-GPP total (RunningMean, τ = `tau_long_term`),
#   `precip_annual` smoothed 1-year precipitation total (RunningMean, τ = `tau_long_term`),
#   `PET_annual`    smoothed 1-year net-radiation potential-evap total (RunningMean,
#                   τ = `tau_long_term`); over `precip_annual` it is the aridity index
#                   driving the online `f0`,
#   `VPDA0_annual`  smoothed 1-year total of VPD·A0 (RunningMean, τ = `tau_long_term`);
#                   over `A0_annual` it is the A0-weighted growing-season VPD,
#   `growing_days`  trailing-year count of growing days, air T > 0 °C (RunningMean,
#                   τ = `tau_long_term`), the online growing-season length,
#   `A0c3_annual`   smoothed 1-year C3 potential-GPP total (RunningMean, τ = `tau_long_term`),
#   `A0c4_annual`   smoothed 1-year C4 potential-GPP total (RunningMean, τ = `tau_long_term`);
#                   the C3/C4 pair drives the online competition for `fractional_c3`,
#   `LAI`           the Zhou et al. (2025) Eq. 16 acclimation lag (RunningMean, τ = 1 day / α),
#                   relaxing toward the instantaneous optimal target `L_opt`.
# The annual totals are a RunningMean of the *annualized* instantaneous rate
# (`year_ref × rate`), so `tau_long_term` sets only the smoothing — the magnitude is a
# one-year total for any `tau_long_term`, as the LAI_max/steady-state formulas require.
# `A0_daily`/`A0_annual`/`precip_annual` set `LAI_max` and the steady-state LAI; the
# cache `area_index.leaf` mirrors the prognostic `LAI` for the rest of the canopy.
function optimal_lai_tivs(component::ZhouOptimalLAIModel{FT}) where {FT}
    seconds_per_day = IP.day(IP.InsolationParameters(FT))
    tau_long_term = component.parameters.tau_long_term
    return (
        ClimaLand.TimeIntegratedVariable(;
            name = :A0_daily,
            reduction = ClimaLand.RunningSum(),
            timescale = seconds_per_day,
        ),
        ClimaLand.TimeIntegratedVariable(;
            name = :A0_annual,
            reduction = ClimaLand.RunningMean(),
            timescale = tau_long_term,
        ),
        ClimaLand.TimeIntegratedVariable(;
            name = :precip_annual,
            reduction = ClimaLand.RunningMean(),
            timescale = tau_long_term,
        ),
        ClimaLand.TimeIntegratedVariable(;
            name = :PET_annual,
            reduction = ClimaLand.RunningMean(),
            timescale = tau_long_term,
        ),
        ClimaLand.TimeIntegratedVariable(;
            name = :VPDA0_annual,
            reduction = ClimaLand.RunningMean(),
            timescale = tau_long_term,
        ),
        ClimaLand.TimeIntegratedVariable(;
            name = :growing_days,
            reduction = ClimaLand.RunningMean(),
            timescale = tau_long_term,
        ),
        ClimaLand.TimeIntegratedVariable(;
            name = :A0c3_annual,
            reduction = ClimaLand.RunningMean(),
            timescale = tau_long_term,
        ),
        ClimaLand.TimeIntegratedVariable(;
            name = :A0c4_annual,
            reduction = ClimaLand.RunningMean(),
            timescale = tau_long_term,
        ),
        ClimaLand.TimeIntegratedVariable(;
            name = :LAI,
            reduction = ClimaLand.RunningMean(),
            timescale = seconds_per_day / component.parameters.alpha,
        ),
    )
end
# `prognostic_vars` (and the matching `_types` / `_domain_names`) are derived from the
# spec metadata in `optimal_lai_tivs`, so the declaration lives in one place.
ClimaLand.prognostic_vars(m::ZhouOptimalLAIModel) =
    ClimaLand.time_integrated_prognostic_vars(optimal_lai_tivs(m))
ClimaLand.prognostic_types(m::ZhouOptimalLAIModel) =
    ClimaLand.time_integrated_prognostic_types(optimal_lai_tivs(m))
ClimaLand.prognostic_domain_names(m::ZhouOptimalLAIModel) =
    ClimaLand.time_integrated_prognostic_domain_names(optimal_lai_tivs(m))

"""
    update_biomass!(
        p,
        Y,
        t,
        component::ZhouOptimalLAIModel{FT},
        canopy,
    ) where {FT}

Updates the optimal-LAI cache: sets SAI and RAI to their prescribed constant
values and mirrors the prognostic `LAI` into the cache area index read by the rest
of the canopy. LAI must be updated here first in `update_aux`, before radiative
transfer reads the area index.
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
    # Climate-responsive f0 from the online aridity index AI = <PET>/<P>, gated by
    # `online_f0` (0 keeps the static artifact f0 already in the cache). Refreshed
    # here in update_aux so compute_LAI_target reads the current f0.
    online = component.parameters.online_f0
    f0_scale = component.parameters.f0_scale
    f0_precip_a = component.parameters.f0_precip_a
    f0_precip_b = component.parameters.f0_precip_b
    f0_precip_c = component.parameters.f0_precip_c
    f0_precip_d = component.parameters.f0_precip_d
    # f0 = f0_source · aridity-targeted f0_scale, where f0_source is the online AI-f0
    # or the raw artifact f0_base. The f0_scale reduction is applied only over the
    # semi-arid precip band [f0_precip_*], leaving the humid tropics/temperate and the
    # deserts unchanged. Recomputed from f0_base each step, so it never compounds;
    # f0_scale=1 is a no-op (== raw f0).
    @. p.canopy.biomass.f0 = aridity_scaled_f0(
        online * f0_from_aridity(
            Y.canopy.biomass.PET_annual,
            Y.canopy.biomass.precip_annual,
        ) + (1 - online) * p.canopy.biomass.f0_base,
        Y.canopy.biomass.precip_annual,
        f0_scale,
        f0_precip_a,
        f0_precip_b,
        f0_precip_c,
        f0_precip_d,
    )
    # When online, recompute vpd_gs as the A0-weighted running mean
    # VPDA0_annual/A0_annual; when static, keep the cache vpd_gs unchanged.
    online_vpd = component.parameters.online_vpd_gs
    @. p.canopy.biomass.vpd_gs =
        online_vpd * (
            Y.canopy.biomass.VPDA0_annual /
            max(Y.canopy.biomass.A0_annual, eps(FT))
        ) + (1 - online_vpd) * p.canopy.biomass.vpd_gs
    # When online, GSL = the trailing-year growing-day count; else keep the
    # cache GSL (static artifact) unchanged.
    online_gsl = component.parameters.online_gsl
    @. p.canopy.biomass.GSL =
        online_gsl * Y.canopy.biomass.growing_days +
        (1 - online_gsl) * p.canopy.biomass.GSL
    # When online, overwrite the P-model's fractional_c3 cache with the dynamic
    # C3/C4 competition value (from the running-mean per-pathway potential GPP);
    # else leave it at the static value the P-model seeded. Read one step later.
    online_c3c4 = component.parameters.online_c3c4
    Mc = canopy.photosynthesis.constants.Mc
    k = component.parameters.k
    @. p.canopy.photosynthesis.fractional_c3 =
        online_c3c4 * c3_fraction_from_competition(
            Y.canopy.biomass.A0c3_annual,
            Y.canopy.biomass.A0c4_annual,
            Mc,
            1 - exp(-k * Y.canopy.biomass.LAI),  # realized fAPAR
        ) + (1 - online_c3c4) * p.canopy.photosynthesis.fractional_c3
    # Mirror the prognostic LAI (a RunningMean relaxing toward the instantaneous
    # optimal target) into the cache area index read by the rest of the canopy.
    @. p.canopy.biomass.area_index.leaf = Y.canopy.biomass.LAI
    # Apply clipping to LAI (same as PrescribedBiomassModel)
    p.canopy.biomass.area_index.leaf .=
        clip.(p.canopy.biomass.area_index.leaf, FT(0.05))
    mask_biomass!(p, Val(canopy.boundary_conditions.prognostic_land_components))
end

"""
    ClimaLand.make_compute_exp_tendency(component::ZhouOptimalLAIModel, canopy)

Advances the optimal-LAI model's four time-integrated variables (declared by
`optimal_lai_tivs`): the 1-day potential-GPP total `A0_daily` as a `RunningSum`, the
smoothed 1-year totals `A0_annual` and `precip_annual` as `RunningMean`s of the
*annualized* instantaneous rate (`year · rate`), and `LAI` as a `RunningMean`
relaxing toward the instantaneous steady-state target `L_opt`:

    dA0_daily/dt      = A0_inst - A0_daily / τ_day,              τ_day  = 1 day,
    dA0_annual/dt     = (year·A0_inst - A0_annual) / τ_long,     τ_long = tau_long_term,
    dprecip_annual/dt = (year·P_inst  - precip_annual) / τ_long,
    dLAI/dt           = (L_opt - LAI) / τ_LAI,                   τ_LAI  = 1 day / α.

`A0_inst` (mol CO2 m^-2 s^-1) is computed once per cache update in `update_biomass!`
and shared by both A0 variables; `P_inst` is the instantaneous precipitation
(mol H2O m^-2 s^-1). Because the annual totals take the mean of the annualized rate,
their steady state is `year · mean(rate)` — a one-year total for any `tau_long_term`,
which sets only the smoothing (a longer value filters the seasonal cycle). Everything
evolves smoothly every timestep with no callback or year-boundary jump.
"""
function ClimaLand.make_compute_exp_tendency(
    component::ZhouOptimalLAIModel{FT},
    canopy,
) where {FT}
    ρ_m_liq = LP.ρ_m_liq(canopy.earth_param_set)  # mol H2O m^-3 (precip volume flux → molar flux)
    # Fixed 365-day reference year (s) that annualizes an instantaneous rate, so the
    # `A0_annual`/`precip_annual` running means hold a one-year total for any `τ`.
    year = FT(365) * IP.day(IP.InsolationParameters(FT))
    specs = optimal_lai_tivs(component)
    # Instantaneous quantity `f` for each variable declared in `optimal_lai_tivs`
    # (a field or `lazy` broadcast). These close over the parent `canopy`, which is
    # only in scope here, so they are supplied to the tendency rather than declared
    # with the variables.
    computes = (;
        A0_daily = (Y, p, t) -> p.canopy.biomass.A0_inst,
        A0_annual = (Y, p, t) -> (@. lazy(year * p.canopy.biomass.A0_inst)),
        # P_liq/P_snow are negative-downward volume fluxes (m/s); negate for a positive total.
        precip_annual = (Y, p, t) -> (@. lazy(
            year * -(p.drivers.P_liq + p.drivers.P_snow) * ρ_m_liq,
        )),
        # Annualized so it matches `precip_annual`'s units and AI = PET/P is dimensionless.
        PET_annual = (Y, p, t) -> begin
            pet = compute_pet_inst(p, canopy)
            @. lazy(year * pet)
        end,
        # The online `vpd_gs` is VPDA0_annual / A0_annual = <VPD·A0>/<A0>, a mean VPD
        # weighted toward the productive period.
        VPDA0_annual = (Y, p, t) -> begin
            vpd_a0 = compute_vpd_a0_inst(p, canopy)
            @. lazy(year * vpd_a0)
        end,
        # Growing-day indicator × 365, so its running mean is the trailing-year count
        # of growing days (air T > 0 C) = the online growing-season length in days.
        growing_days = (Y, p, t) ->
            (@. lazy(ifelse(p.drivers.T > FT(273.15), FT(365), FT(0)))),
        # Pure-C3 and pure-C4 potential GPP, annualized; the pair drives the C3/C4
        # competition through the advantage (A0c4 - A0c3)/A0c3.
        A0c3_annual = (Y, p, t) -> begin
            a0_c3 = compute_A0_inst(p, canopy, one(FT))
            @. lazy(year * a0_c3)
        end,
        A0c4_annual = (Y, p, t) -> begin
            a0_c4 = compute_A0_inst(p, canopy, zero(FT))
            @. lazy(year * a0_c4)
        end,
        LAI = (Y, p, t) -> compute_LAI_target(Y, p, canopy),
    )
    function compute_exp_tendency!(dY, Y, p, t)
        # Compute the instantaneous potential GPP once (drivers, `par_d`, `βm` fresh
        # after `update_aux`); both A0 reductions read it.
        compute_A0_inst!(p.canopy.biomass.A0_inst, p, canopy)
        ClimaLand.time_integrated_tendency!(
            dY.canopy.biomass,
            Y.canopy.biomass,
            Y,
            p,
            t,
            specs,
            computes,
        )
    end
    return compute_exp_tendency!
end

"""
    set_historical_cache!(p, Y0, model::ZhouOptimalLAIModel, canopy)

Copies the optimal-LAI model's static spatially-varying inputs into the cache, where
the `LAI_max` and steady-state-LAI formulas read them: the growing season length
`GSL` (days), the growing-season mean VPD `vpd_gs` (Pa) used in the water-limitation
WUE factor, and the fraction of precipitation available for transpiration `f0`. They
are taken from `model.optimal_lai_inputs`, typically built by
`optimal_lai_initial_conditions`, and are rebuilt at every setup, including on
restart.

The model's prognostic state (`LAI`, `A0_daily`, `A0_annual`, `precip_annual`) is not
touched here; it is initialized with the rest of `Y` by
`ClimaLand.Simulations.set_canopy_biomass_initial_conditions!` in `set_ic!`.
"""
function set_historical_cache!(p, Y0, model::ZhouOptimalLAIModel, canopy)
    optimal_lai_inputs = model.optimal_lai_inputs
    # `f0_base` holds the raw artifact f0; `update_biomass!` recomputes the f0 actually
    # used as f0_source scaled over the semi-arid band each step, so nothing compounds.
    p.canopy.biomass.GSL .= optimal_lai_inputs.GSL
    p.canopy.biomass.vpd_gs .= optimal_lai_inputs.vpd_gs
    p.canopy.biomass.f0_base .= optimal_lai_inputs.f0
    p.canopy.biomass.f0 .= optimal_lai_inputs.f0
    return nothing
end

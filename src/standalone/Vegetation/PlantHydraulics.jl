module PlantHydraulics
using ClimaLand
using ClimaLand.Soil
using ClimaUtilities.TimeVaryingInputs
import ClimaUtilities.TimeVaryingInputs:
    TimeVaryingInput, AbstractTimeVaryingInput
import ClimaUtilities.TimeManager: ITime
import NCDatasets, ClimaCore, Interpolations # Needed to load TimeVaryingInputs
using ..ClimaLand.Canopy: AbstractCanopyComponent
using ClimaLand: PrescribedGroundConditions
using ClimaCore
import ClimaLand.Parameters as LP
import ClimaParams as CP
using DocStringExtensions
using LazyBroadcast: lazy

import ClimaLand:
    make_update_aux,
    make_compute_exp_tendency,
    prognostic_vars,
    prognostic_types,
    auxiliary_vars,
    auxiliary_domain_names,
    prognostic_domain_names,
    name,
    total_liq_water_vol_per_area!
export PlantHydraulicsModel,
    AbstractPlantHydraulicsModel,
    water_flux,
    effective_saturation,
    augmented_liquid_fraction,
    water_retention_curve,
    inverse_water_retention_curve,
    root_water_flux_per_ground_area!,
    PlantHydraulicsParameters,
    AbstractConductivityModel,
    AbstractRetentionModel,
    LinearRetentionCurve,
    Weibull,
    hydraulic_conductivity

"""
    AbstractPlantHydraulicsModel{FT} <: AbstractCanopyComponent{FT}

An abstract type for plant hydraulics models.
"""
abstract type AbstractPlantHydraulicsModel{FT} <: AbstractCanopyComponent{FT} end

ClimaLand.name(::AbstractPlantHydraulicsModel) = :hydraulics

"""
    PlantHydraulicsParameters

A struct for holding parameters of the PlantHydraulics Model.
$(DocStringExtensions.FIELDS)
"""
struct PlantHydraulicsParameters{FT <: AbstractFloat, CP, RP}
    "porosity (m3/m3)"
    ν::FT
    "storativity (m3/m3)"
    S_s::FT
    "Conductivity model and parameters"
    conductivity_model::CP
    "Water retention model and parameters"
    retention_model::RP
end

"""
    PlantHydraulicsParameters(;
        ν::FT,
        S_s::FT,
        conductivity_model,
        retention_model,
    )

Constructor for PlantHydraulicsParameters.
"""
function PlantHydraulicsParameters(;
    ν::FT,
    S_s::FT,
    conductivity_model,
    retention_model,
) where {FT}
    return PlantHydraulicsParameters{
        FT,
        typeof(conductivity_model),
        typeof(retention_model),
    }(
        ν,
        S_s,
        conductivity_model,
        retention_model,
    )
end

function PlantHydraulicsParameters(
    toml_dict::CP.ParamDict;
    ν = toml_dict["plant_nu"],
    S_s = toml_dict["plant_S_s"],
    conductivity_model,
    retention_model,
)
    FT = typeof(ν)
    return PlantHydraulicsParameters{
        FT,
        typeof(conductivity_model),
        typeof(retention_model),
    }(
        ν,
        S_s,
        conductivity_model,
        retention_model,
    )
end


"""
    PlantHydraulicsModel{FT, PS, T, AA} <: AbstractPlantHydraulicsModel{FT}

Defines, and constructs instances of, the PlantHydraulicsModel type, which is used
for simulation flux of water to/from soil, along roots of different depths,
along a stem, to a leaf, and ultimately being lost from the system by
transpiration. Note that the canopy height is specified as part of the
PlantHydraulicsModel and the biomass model, these must be consistent.

The model can be used in Canopy standalone mode by prescribing
the soil matric potential at the root tips or flux in the roots. There is also the
option (intendend only for debugging) to use a prescribed transpiration rate.

$(DocStringExtensions.FIELDS)
"""
struct PlantHydraulicsModel{FT, PS, AA <: AbstractArray{FT}} <:
       AbstractPlantHydraulicsModel{FT}
    "The number of stem compartments for the plant; can be zero"
    n_stem::Int64
    "The number of leaf compartments for the plant; must be >=1"
    n_leaf::Int64
    "The height of the center of each leaf compartment/stem compartment, in meters"
    compartment_midpoints::AA
    "The height of the compartments' top faces, in meters. The canopy height is the last element of the vector."
    compartment_surfaces::AA
    "The label (:stem or :leaf) of each compartment"
    compartment_labels::Vector{Symbol}
    "Parameters required by the Plant Hydraulics model"
    parameters::PS
end

function PlantHydraulicsModel{FT}(;
    n_stem::Int64,
    n_leaf::Int64,
    compartment_midpoints::Vector{FT},
    compartment_surfaces::Vector{FT},
    parameters::PlantHydraulicsParameters{FT},
) where {FT}
    @assert (n_leaf + n_stem) == length(compartment_midpoints)
    @assert (n_leaf + n_stem) + 1 == length(compartment_surfaces)
    for i in 1:length(compartment_midpoints)
        @assert compartment_midpoints[i] ==
                ((compartment_surfaces[i + 1] - compartment_surfaces[i]) / 2) +
                compartment_surfaces[i]
    end
    compartment_labels = Vector{Symbol}(undef, n_stem + n_leaf)
    for i in 1:(n_stem + n_leaf)
        if i <= n_stem
            compartment_labels[i] = :stem
        else
            compartment_labels[i] = :leaf
        end
    end
    return PlantHydraulicsModel{
        FT,
        typeof(parameters),
        typeof(compartment_midpoints),
    }(
        n_stem,
        n_leaf,
        compartment_midpoints,
        compartment_surfaces,
        compartment_labels,
        parameters,
    )
end

"""
    prognostic_vars(model::PlantHydraulicsModel)

A function which returns the names of the prognostic
variables of the `PlantHydraulicsModel`.
"""
prognostic_vars(model::PlantHydraulicsModel) = (:ϑ_l,)

"""
    auxiliary_vars(model::PlantHydraulicsModel)

A function which returns the names of the auxiliary
variables of the `PlantHydraulicsModel`,
the water potential `ψ` (m), the volume flux\\*\\cross section `fa` (1/s),
and the volume flux\\*root cross section in the roots `fa_roots` (1/s),
where the cross section can be represented by an area index.

The water potential ψ is of length `n`, where `n` is the number of compartments,
and `fa` is of length `n-1`. If `n=1`, `fa` is not included. `fa_roots` is a scalar,
the bottom boundary flux.
"""
function auxiliary_vars(model::PlantHydraulicsModel)
    n = model.n_stem + model.n_leaf
    if n > 1
        return (:ψ, :fa, :fa_roots)
    else
        return (:ψ, :fa_roots)
    end
end


"""
    ClimaLand.prognostic_types(model::PlantHydraulicsModel{FT}) where {FT}

Defines the prognostic types for the PlantHydraulicsModel.
"""
ClimaLand.prognostic_types(model::PlantHydraulicsModel{FT}) where {FT} =
    (NTuple{model.n_stem + model.n_leaf, FT},)
ClimaLand.prognostic_domain_names(::PlantHydraulicsModel) = (:surface,)

"""
    ClimaLand.auxiliary_types(model::PlantHydraulicsModel{FT}) where {FT}

Defines the auxiliary types for the PlantHydraulicsModel.
"""
function ClimaLand.auxiliary_types(model::PlantHydraulicsModel{FT}) where {FT}
    n = model.n_stem + model.n_leaf
    if n > 1
        return (
            NTuple{model.n_stem + model.n_leaf, FT},
            NTuple{model.n_stem + model.n_leaf - 1, FT},
            FT,
        )
    else
        return (NTuple{model.n_stem + model.n_leaf, FT}, FT)
    end
end

function ClimaLand.auxiliary_domain_names(model::PlantHydraulicsModel)
    n = model.n_stem + model.n_leaf
    if n > 1
        return (:surface, :surface, :surface)
    else
        return (:surface, :surface)
    end
end

"""
    harmonic_mean(x::FT,y::FT) where {FT}

Computes the harmonic mean of x >=0 and y >=0; returns zero if either
x or y are zero.
"""
harmonic_mean(x::FT, y::FT) where {FT} = x * y / max(x + y, eps(FT))

"""
    water_flux(
        z1,
        z2,
        ψ1,
        ψ2,
        K1,
        K2,
    ) where {FT}

Computes the water flux given the absolute potential ψ (pressure/(ρg))
 and the conductivity K (m/s) at the center of the two layers
with midpoints z1 and z2.

We currently assuming a harmonic
mean for effective conducticity between the two layers
(see CLM Technical Documentation).

To account for different path lengths in the two compartments Δz1 and
Δz2, we would require the following conductance k (1/s)
k\\_eff = K1/Δz1*K2/Δz2/(K1/Δz1+K2/Δz2)
and a water flux of
F = -k\\_eff * (ψ1 +z1 - ψ2 - z2) (m/s).

This currently assumes the path lengths are equal.
"""
function water_flux(z1::FT, z2::FT, ψ1::FT, ψ2::FT, K1::FT, K2::FT) where {FT}
    K_eff = harmonic_mean(K1, K2)
    flux = -K_eff * ((ψ2 - ψ1) / (z2 - z1) + 1)
    return flux # (m/s)
end

"""
    AbstractConductivityModel{FT <: AbstractFloat}

An abstract type for the plant hydraulics conductivity model.
"""
abstract type AbstractConductivityModel{FT <: AbstractFloat} end

Base.broadcastable(x::AbstractConductivityModel) = tuple(x)
"""
    AbstractRetentionModel{FT <: AbstractFloat}

An abstract type for the plant retention curve model.
"""
abstract type AbstractRetentionModel{FT <: AbstractFloat} end

Base.broadcastable(x::AbstractRetentionModel) = tuple(x)


"""
    Weibull{FT} <: AbstractConductivityModel{FT}

A concrete type specifying that a Weibull conductivity model is to be used;
the struct contains the require parameters for this model.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Weibull{FT} <: AbstractConductivityModel{FT}
    "Maximum Water conductivity in the above-ground plant compartments (m/s) at saturation"
    K_sat::FT
    "The absolute water potential in xylem (or xylem water potential) at which ∼63%
    of maximum xylem conductance is lost (Liu, 2020)."
    ψ63::FT
    "Weibull parameter c, which controls shape the shape of the conductance curve (Sperry, 2016)."
    c::FT
end

function Weibull(
    toml_dict::CP.ParamDict;
    K_sat_plant = toml_dict["K_sat_plant"],
    ψ63 = toml_dict["psi_63"],
    c = toml_dict["Weibull_c"],
)
    FT = CP.float_type(toml_dict)
    return Weibull{FT}(K_sat_plant, ψ63, c)
end

"""
    hydraulic_conductivity(conductivity_params::Weibull{FT}, ψ::FT) where {FT}

Computes the hydraulic conductivity at a point, using the
Weibull formulation, given the potential ψ.
"""
function hydraulic_conductivity(
    conductivity_params::Weibull{FT},
    ψ::FT,
) where {FT}
    (; K_sat, ψ63, c) = conductivity_params
    if ψ <= FT(0)
        K = exp(-(ψ / ψ63)^c)
    else
        K = FT(1)
    end
    return K * K_sat # (m/s)
end

"""
    LinearRetentionCurve{FT} <: AbstractRetentionModel{FT}

A concrete type specifying that a linear water retention  model is to be used;
the struct contains the require parameters for this model.

When ψ = 0, the effective saturation is one, so the intercept
is not a free parameter, and only the slope must be specified.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct LinearRetentionCurve{FT} <: AbstractRetentionModel{FT}
    "Bulk modulus of elasticity and slope of potential to volume curve. See also Corcuera, 2002, and Christoffersen, 2016."
    a::FT
end

function LinearRetentionCurve(toml_dict::CP.ParamDict; a = toml_dict["a"])
    FT = CP.float_type(toml_dict)
    return LinearRetentionCurve{FT}(a)
end

"""
    water_retention_curve(
        S_l::FT,
        b::FT,
        ν::FT,
        S_s::FT) where {FT}

Returns the potential ψ given the effective saturation S at a point, according
to a linear model for the retention curve with parameters specified
by `retention_params`.
"""
function water_retention_curve(
    retention_params::LinearRetentionCurve{FT},
    S_l::FT,
    ν::FT,
    S_s::FT,
) where {FT}
    (; a) = retention_params
    if S_l <= FT(1)
        ψ = 1 / a * (S_l - 1) # ψ(S_l=1)=0.
    else
        ϑ_l = augmented_liquid_fraction(ν, S_l)
        ψ = (ϑ_l - ν) / S_s
    end
    return ψ # (m)
end

"""
    inverse_water_retention_curve(
        ψ::FT,
        b::FT,
        ν::FT,
        S_s::FT) where {FT}

Returns the effective saturation given the potential at a point, according
to the linear retention curve model.
"""
function inverse_water_retention_curve(
    retention_params::LinearRetentionCurve{FT},
    ψ::FT,
    ν::FT,
    S_s::FT,
) where {FT}
    (; a) = retention_params
    if ψ <= FT(0)
        S_l = ψ * a + 1 # ψ(S_l=1)=0.
    else
        ϑ_l = ψ * S_s + ν
        S_l = effective_saturation(ν, ϑ_l)
    end
    return S_l #(m3/m3)
end

"""
    augmented_liquid_fraction(
        ν::FT,
        S_l::FT) where {FT}

Computes the augmented liquid fraction from porosity and
effective saturation.

Augmented liquid fraction allows for
oversaturation: an expansion of the volume of space
available for storage in a plant compartment.
"""
function augmented_liquid_fraction(ν::FT, S_l::FT) where {FT}
    ϑ_l = S_l * ν # ϑ_l can be > ν
    safe_ϑ_l = max(ϑ_l, eps(FT))
    return safe_ϑ_l # (m3 m-3)
end

"""
    effective_saturation(
        ν::FT,
        ϑ_l::FT) where {FT}

Computes the effective saturation given the augmented liquid fraction.
"""
function effective_saturation(ν::FT, ϑ_l::FT) where {FT}
    S_l = ϑ_l / ν # S_l can be > 1
    safe_S_l = max(S_l, eps(FT))
    return safe_S_l # (m3 m-3)
end

"""
    make_compute_exp_tendency(model::PlantHydraulicsModel, _)

A function which creates the compute_exp_tendency! function for the PlantHydraulicsModel.
The compute_exp_tendency! function must comply with a rhs function of SciMLBase.jl.

Below, `fa` denotes a flux multiplied by the relevant cross section
(per unit area ground, or area index, AI). The tendency for the
ith compartment can be written then as:
∂ϑ[i]/∂t = 1/(AI*dz)[fa[i]-fa[i+1]).

Note that if the area_index is zero because no plant is present,
AIdz is zero, and the fluxes `fa` appearing in the numerator are
zero because they are scaled by AI.

To prevent dividing by zero, we change AI/(AI x dz)" to
"AI/max(AI x dz, eps(FT))"
"""
function make_compute_exp_tendency(
    model::PlantHydraulicsModel{FT},
    canopy,
) where {FT}
    function compute_exp_tendency!(dY, Y, p, t)
        area_index = p.canopy.biomass.area_index
        n_stem = model.n_stem
        n_leaf = model.n_leaf
        fa_roots = p.canopy.hydraulics.fa_roots

        # Inside of a loop, we need to use a single dollar sign
        # for indexing into Fields of Tuples in non broadcasted
        # expressions, and two dollar signs for
        # for broadcasted expressions using the macro @.
        # field.:($index) .= value # works
        # @ field.:($$index) = value # works
        i_end = n_stem + n_leaf
        if i_end == 1 # single compartment
            AI = getproperty(area_index, model.compartment_labels[1]) # this is a field; should not allocate here
            dz = model.compartment_surfaces[2] - model.compartment_surfaces[1] # currently this is a scalar. in the future, this will be a field.
            @. dY.canopy.hydraulics.ϑ_l.:1 =
                1 / max(AI * dz, eps(FT)) *
                (fa_roots - p.canopy.turbulent_fluxes.vapor_flux)
        else # multiple layers
            fa = p.canopy.hydraulics.fa
            @inbounds for i in 1:(n_stem + n_leaf)
                im1 = i - 1
                ip1 = i + 1
                # To prevent dividing by zero, change AI/(AI x dz)" to AI/max(AI x dz, eps(FT))"
                AI = getproperty(area_index, model.compartment_labels[i]) # this is a field; should not allocate here
                dz =
                    model.compartment_surfaces[ip1] -
                    model.compartment_surfaces[i] # currently this is a scalar. in the future, this will be a field.
                if i == 1
                    @inbounds @. dY.canopy.hydraulics.ϑ_l.:($$i) =
                        1 / max(AI * dz, eps(FT)) * (fa_roots - fa.:($$i))
                elseif i == i_end
                    @inbounds @. dY.canopy.hydraulics.ϑ_l.:($$i) =
                        1 / max(AI * dz, eps(FT)) *
                        (fa.:($$im1) - p.canopy.turbulent_fluxes.vapor_flux)
                else
                    @inbounds @. dY.canopy.hydraulics.ϑ_l.:($$i) =
                        1 / max(AI * dz, eps(FT)) * (fa.:($$im1) - fa.:($$i))
                end
            end
        end
    end
    return compute_exp_tendency!
end

"""
    root_water_flux_per_ground_area!(
        fa::ClimaCore.Fields.Field,
        ground::PrescribedGroundConditions,
        model::PlantHydraulicsModel{FT},
        canopy,
        Y::ClimaCore.Fields.FieldVector,
        p::NamedTuple,
        t,
    ) where {FT}

A method which computes the water flux between the soil and the stem, via the roots,
and multiplied by the RAI, in the case of a model running without a prognostic
soil model:

Flux = -K_eff x [(ψ_stem - ψ_soil)/(z_stem - z_soil) + 1], where
K_eff = K_soil K_stem /(K_stem + K_soil)

Note that in `PrescribedSoil` mode, we compute the flux using K_soil = K_plant(ψ_soil)
and K_stem = K_plant(ψ_stem). In `PrognosticSoil` mode, we compute the flux using
K_soil = K_soil(ψ_soil) and K_stem = K_plant(ψ_stem). The latter is a better model, but
our `PrescribedSoil` struct does not store K_soil, only θ_soil.
θ_soil is converted to ψ_soil using the retention curve supplied by hydrology_cm.

The returned flux is per unit ground area. This assumes that the stem compartment
is the first element of `Y.canopy.hydraulics.ϑ_l`.
"""
function root_water_flux_per_ground_area!(
    fa::ClimaCore.Fields.Field,
    ground::PrescribedGroundConditions,
    model::PlantHydraulicsModel{FT},
    canopy,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    t,
) where {FT}
    rooting_depth = canopy.biomass.rooting_depth
    (; conductivity_model,) = model.parameters
    area_index = p.canopy.biomass.area_index
    # We can index into a field of Tuple{FT} to extract a field of FT
    # using the following notation: field.:index
    ψ_base = p.canopy.hydraulics.ψ.:1

    # compute the soil water potential from soil water content from retention curve
    soil_saturation = @. lazy(
        max(
            min(
                Soil.effective_saturation(ground.ν, p.drivers.θ, ground.θ_r),
                1,
            ),
            eps(FT),
        ),
    )
    ψ_soil = @. lazy(matric_potential(ground.hydrology_cm, soil_saturation))

    above_ground_area_index =
        harmonic_mean.(
            getproperty(area_index, model.compartment_labels[1]),
            getproperty(area_index, :root),
        )
    # since rooting_depth is positive by convention, add the sign in here to
    # convert it to a coordinate: z_roots = -rooting_depth
    @. fa .=
        water_flux(
            -rooting_depth,
            model.compartment_midpoints[1],
            ψ_soil,
            ψ_base,
            hydraulic_conductivity(conductivity_model, ψ_soil),
            hydraulic_conductivity(conductivity_model, ψ_base),
        ) * above_ground_area_index
end

"""
    ClimaLand.total_liq_water_vol_per_area!(
        surface_field,
        model::PlantHydraulicsModel,
        Y,
        p,
        t,
)

A function which updates `surface_field` in place with the value of
the plant hydraulic models total water volume.

Note that this is general for any number of canopy layers, but it assumes
that the LAI and SAI given are per layer. This is distinct from the BigLeaf
approach, in which the LAI and SAI refer to the integrated area index with heigh.
"""
function ClimaLand.total_liq_water_vol_per_area!(
    surface_field,
    model::PlantHydraulicsModel,
    Y,
    p,
    t,
)
    area_index = p.canopy.biomass.area_index
    dz =
        model.compartment_surfaces[2:end] .-
        model.compartment_surfaces[1:(end - 1)]
    labels = model.compartment_labels
    surface_field .= 0
    n = length(labels)
    for i in 1:n
        surface_field .+=
            dz[i] .* getproperty(area_index, labels[i]) .*
            Y.canopy.hydraulics.ϑ_l.:($i)
    end
    return nothing
end


"""
   update_hydraulics!(p, Y, hydraulics::PlantHydraulicsModel, canopy)

Updates the following cache variables in place:
- p.canopy.hydraulics.ψ
- p.canopy.hydraulics.fa[1:end-1] (within plant fluxes)

Other types of AbstractPlantHydraulicsModel may update different variables.
"""
function update_hydraulics!(p, Y, hydraulics::PlantHydraulicsModel, canopy)
    ψ = p.canopy.hydraulics.ψ
    ϑ_l = Y.canopy.hydraulics.ϑ_l
    area_index = p.canopy.biomass.area_index
    n_stem = hydraulics.n_stem
    n_leaf = hydraulics.n_leaf
    (; retention_model, conductivity_model, S_s, ν) = hydraulics.parameters
    # We can index into a field of Tuple{FT} to extract a field of FT
    # using the following notation: field.:index
    @inbounds @. ψ.:1 = PlantHydraulics.water_retention_curve(
        retention_model,
        PlantHydraulics.effective_saturation(ν, ϑ_l.:1),
        ν,
        S_s,
    )
    # Inside of a loop, we need to use a single dollar sign
    # for indexing into Fields of Tuples in non broadcasted
    # expressions, and two dollar signs for
    # for broadcasted expressions using the macro @.
    # field.:($index) .= value # works
    # @ field.:($$index) = value # works
    @inbounds for i in 1:(n_stem + n_leaf - 1)
        ip1 = i + 1
        @. ψ.:($$ip1) = PlantHydraulics.water_retention_curve(
            retention_model,
            PlantHydraulics.effective_saturation(ν, ϑ_l.:($$ip1)),
            ν,
            S_s,
        )

        areai = getproperty(area_index, hydraulics.compartment_labels[i])
        areaip1 = getproperty(area_index, hydraulics.compartment_labels[ip1])

        # Compute the flux*area between the current compartment `i`
        # and the compartment above.
        @. p.canopy.hydraulics.fa.:($$i) =
            PlantHydraulics.water_flux(
                hydraulics.compartment_midpoints[i],
                hydraulics.compartment_midpoints[ip1],
                ψ.:($$i),
                ψ.:($$ip1),
                PlantHydraulics.hydraulic_conductivity(
                    conductivity_model,
                    ψ.:($$i),
                ),
                PlantHydraulics.hydraulic_conductivity(
                    conductivity_model,
                    ψ.:($$ip1),
                ),
            ) * PlantHydraulics.harmonic_mean(areaip1, areai)
    end
end
end

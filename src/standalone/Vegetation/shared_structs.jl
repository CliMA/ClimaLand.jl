export SpectralDiscretization

using StaticArrays

"""
    SpectralDiscretization

A type used to represent the spectral discretization of radiation models.
Consists of the wavelength boundaries between each band in a discretization.
"""
struct SpectralDiscretization{FT <: AbstractFloat}
    "Wavelength boundaries between each band in the discretization (m)"
    λ::Tuple{Vararg{FT}}
    "Banded irradiance curve for the solar radiation (W/m^2)"
    I::Tuple{Vararg{FT}}
    "Proportions of PAR radiation in each band"
    PAR_proportions::Tuple{Vararg{FT}}
    "Proportions of NIR radiation in each band"
    NIR_proportions::Tuple{Vararg{FT}}
end

"""
    SolarIrradianceCurve(λ::FT) where {FT}

Function to approximate the solar irradiance curve at a given wavelength. The
solar irradiance curve gives the intensity of solar radiation over the spectrum.
"""
function SolarIrradianceCurve(λ::FT) where {FT}
    return (2 * 6.26e-34 * 3e8 / λ^5) * (1 / (exp(6.26e-34 * 3e8 / (λ * 1.381e-23 * 5778)) - 1))
end

"""
    SolarIrradianceArea(λ1::FT, λ2::FT) where{FT}

Quickly computes the area under the solar irradiance curve between two
wavelengths using a step size of 1 nm.
"""
function SolarIrradianceArea(λ1::FT, λ2::FT) where {FT}
    # Create an mvector of the steps
    num_steps = Int(round((λ2 - λ1) * 1e9))
    # Compute the area via the trapezoid method
    return sum(SolarIrradianceCurve.(range(λ1, stop=λ2, length=num_steps))) * (λ2 - λ1) / num_steps
end

"""
    ComputeBandedIrradianceCurve(
        λ::Tuple{Vararg{FT}},
    ) where {FT}

Assigns a proportion of total incoming solar radiation to each band in the
spectral discretization to allow the computation of the input radiation for
a discretization with an arbitrary number of bands.
"""
function ComputeBandedIrradianceCurve(λ::Tuple{Vararg{FT}}) where {FT}
    # Use the trapezoid method to compute the area under the curve between each
    trap_sums = MVector{length(λ) - 1, FT}(zeros(length(λ) - 1)...)
    for i in 1:length(λ) - 1
        trap_sums[i] = SolarIrradianceArea(λ[i], λ[i + 1])
    end
    # Now take each band's area as a proportion of the total area under the
    # curve
    return Tuple(trap_sums ./ sum(trap_sums))
end

"""
    ComputeBandedProportions(
        λ::Tuple{Vararg{FT}},
    ) where {FT}

Use the solar irradiance curve to compute the proportion of each band in a
spectral discretization that is within a specific range, used for PAR and NIR
radiation.
"""
function ComputeBandedProportions(λ::Tuple{Vararg{FT}}, bounds::Tuple{FT, FT}) where {FT}
    # PAR proportions for each band
    proportions = MVector{length(λ) - 1, FT}(zeros(length(λ) - 1)...)

    # In each band, determine whether the band is entirely within the PAR range,
    # partly in the PAR range, or entirely outside the PAR range
    for i in 1:length(λ) - 1
        # If the band is entirely within the range, the proportion of radiation
        # in the range is 1
        if λ[i] >= bounds[1] && λ[i + 1] <= bounds[2]
            proportions[i] = 1
        # If the band is entirely outside the range, the proportion of radiation
        # in the range is 0
        elseif λ[i] >= bounds[2] || λ[i + 1] <= bounds[1]
            proportions[i] = 0
        # If the band is partly in the range, we need to know what part of 
        # the band is in the range
        elseif λ[i] <= bounds[1] && λ[i + 1] <= bounds[2]
            # trapezoidal area of the band that is in the PAR range
            bound_trap = SolarIrradianceArea(bounds[1], λ[i + 1])
            # trapezoidal area of the band
            trap = SolarIrradianceArea(λ[i], λ[i + 1])
            # proportion of PAR radiation in the band
            proportions[i] = bound_trap / trap
        else
            # trapezoidal area of the band that is in the PAR range
            bound_trap = SolarIrradianceArea(λ[i], bounds[2])
            println(bound_trap)
            # trapezoidal area of the band
            trap = SolarIrradianceArea(λ[i], λ[i + 1])
            # proportion of PAR radiation in the band
            proportions[i] = bound_trap / trap
        end
    end
    return Tuple(proportions)
end

"""
    ComputeBandedPAR(
        λ::Tuple{Vararg{FT}},
    ) where {FT}

Use the solar irradiance curve to compute the proportion of each band in a
spectral discretization that is within the PAR range.
"""
function ComputeBandedPAR(λ::Tuple{Vararg{FT}}) where {FT}
    return ComputeBandedProportions(λ, (400e-9, 700e-9))
end

"""
    ComputeBandedNIR(
        λ::Tuple{Vararg{FT}},
    ) where {FT}
Use the solar irradiance curve to compute the proportion of each band in a
spectral discretization that is within the NIR range.
"""
function ComputeBandedNIR(λ::Tuple{Vararg{FT}}) where {FT}
    return ComputeBandedProportions(λ, (700e-9, 2500e-9))
end

"""
    SpectralDiscretization(
        λ::Tuple{Vararg{FT}},
    ) where {FT}
    
Creates a spectral discretization with the given wavelength boundaries.
"""
function SpectralDiscretization(λ::Tuple{Vararg{FT}}) where {FT}
    # Check that the given boundaries are in ascending order and cover the full
    # SW range from 100 nm to 3000 nm
    @assert λ[1] == 100e-9 && λ[end] == 3000e-9
    @assert all(λ[i] < λ[i + 1] for i in 1:length(λ) - 1)
    I = ComputeBandedIrradianceCurve(λ)
    PAR_proportions = ComputeBandedPAR(λ)
    NIR_proportions = ComputeBandedNIR(λ)
    return SpectralDiscretization{FT}(λ, I, PAR_proportions, NIR_proportions)
end
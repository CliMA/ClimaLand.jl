"""
    HeatParameters{FT <: AbstractFloat}

A parameter structure for the heat equation. In this simplest
form, we assume the conductivity and volumetric heat capacity
of the soil are constant.
$(DocStringExtensions.FIELDS)
"""
struct HeatParameters{FT <: AbstractFloat}
    "The thermal conductivity, W/m/K"
    κ::FT
    "The volumetric heat capacity, J/m^3/K"
    ρc_s::FT
end

"""
    SoilEnergyHydrology <: AbstractSoilModel

A model for simulating the flow of water and heat 
in a porous medium by solving the Richardson-Richards equation
and the heat equation, including terms for phase change.

$(DocStringExtensions.FIELDS)
"""
struct SoilEnergyHydrology{FT, PSR, PSH, PS, D, C, BCR, BCH, S} <:
       AbstractSoilModel{FT}
    "the parameter set for RR equation"
    hydrology_param_set::PSR
    "the heat equation parameter set"
    thermal_param_set::PSH
    "Physical constants"
    earth_param_set::PS
    "the soil domain, using ClimaCore.Domains"
    domain::D
    "the domain coordinates"
    coordinates::C
    "the boundary conditions for RRE, of type AbstractSoilBoundaryConditions"
    hydrology_boundary_conditions::BCR
    "the boundary conditions for heat equation, of type AbstractSoilBoundaryConditions"
    thermal_boundary_conditions::BCH
    "A tuple of sources, each of type AbstractSoilSource"
    sources::S
end

"""
    SoilEnergyHydrology{FT}(;
        hydrology_param_set::PSR,
        thermal_param_set::PSH,
        earth_param_set::AbstractEarthParameterSet,
        domain::D,
        hydrology_boundary_conditions::AbstractSoilBoundaryConditions{FT},
        thermal_boundary_conditions::AbstractSoilBoundaryConditions{FT},
        sources::Tuple,
    ) where {FT, PSR, PSH, D}

A constructor for a `SoilEnergyHydrology` model.
"""
function SoilEnergyHydrology{FT}(;
                                 hydrology_param_set::RichardsParameters{FT},
                                 thermal_param_set::HeatParameters{FT},
                                 earth_param_set::AbstractEarthParameterSet,
                                 domain::D,
                                 hydrology_boundary_conditions::AbstractSoilBoundaryConditions{FT},
                                 thermal_boundary_conditions::AbstractSoilBoundaryConditions{FT},
                                 sources::Tuple,
                                 ) where {FT, D}
    coords = coordinates(domain)
    args = (
        hydrology_param_set,
        thermal_param_set,
        earth_param_set,
        domain,
        coords,
        hydrology_boundary_conditions,
        thermal_boundary_conditions,
        sources,
    )
    SoilEnergyHydrology{FT, typeof.(args)...}(args...)
end

"""
    make_rhs(model::SoilEnergyHydrology)

An extension of the function `make_rhs`, for the integrated soil
energy and heat equations, including phase change.

This function creates and returns a function which computes the entire
right hand side of the PDE for `Y.soil.ϑ_l, Y.soil.θ_i, Y.soil.ρe_int`, 
and updates `dY.soil` in place with those values.

This has been written so as to work with Differential Equations.jl.
"""
function ClimaLSM.make_rhs(model::SoilEnergyHydrology{FT}) where {FT}
    function rhs!(dY, Y, p, t)
        @unpack ν, vg_α, vg_n, vg_m, Ksat, S_s, θ_r = model.hydrology_param_set
        @unpack κ = model.thermal_param_set
        ρe_int_l = volumetric_internal_energy_liq.(p.soil.T, Ref(model.earth_param_set))
        z = model.coordinates.z


        hydrology_top_flux_bc, hydrology_bot_flux_bc =
            boundary_fluxes(model.hydrology_boundary_conditions, p, t)
        thermal_top_flux_bc, thermal_bot_flux_bc =
            boundary_fluxes(model.thermal_boundary_conditions, p, t)

        interpc2f = Operators.InterpolateC2F()
        gradc2f = Operators.GradientC2F()
        # Richards-Richardson RHS:

        divf2c_rre = Operators.DivergenceF2C(
            top = Operators.SetValue(Geometry.WVector(hydrology_top_flux_bc)),
            bottom = Operators.SetValue(Geometry.WVector(hydrology_bot_flux_bc)),
        )
        @. dY.soil.ϑ_l =
            -(divf2c_rre(-interpc2f(p.soil.K) * gradc2f(p.soil.ψ + z)))
        dY.soil.θ_i .= ClimaCore.Fields.zeros(FT, axes(Y.soil.θ_i))

        # Heat equation RHS
        divf2c_heat = Operators.DivergenceF2C(
            top = Operators.SetValue(Geometry.WVector(thermal_top_flux_bc)),
            bottom = Operators.SetValue(Geometry.WVector(thermal_bot_flux_bc)),
        )
        @. dY.soil.ρe_int =
            -divf2c_heat(
                -κ * gradc2f(p.soil.T) -
                interpc2f(ρe_int_l * p.soil.K) * gradc2f(p.soil.ψ + z),
            )
        # Horizontal contributions
        horizontal_components!(dY, model.domain, model, p)

        for src in model.sources
            source!(dY, src, Y, p)
        end

        # This has to come last
        dss!(dY, model.domain)
    end
    return rhs!
end

"""
   horizontal_components!(dY::ClimaCore.Fields.FieldVector,
                          domain::HybridBox,
                          model::SoilEnergyHydrology,
                          p::ClimaCore.Fields.FieldVector)

Updates dY in place by adding in the tendency terms resulting from
horizontal derivative operators.

In the case of a hybrid box domain, the horizontal contributions are
computed using the WeakDivergence and Gradient operators.
"""
function horizontal_components!(
    dY::ClimaCore.Fields.FieldVector,
    domain::HybridBox,
    model::SoilEnergyHydrology,
    p::ClimaCore.Fields.FieldVector,
)
    κ = model.thermal_param_set.κ
    ρe_int_l = volumetric_internal_energy_liq.(p.soil.T, Ref(model.earth_param_set))

    hdiv = Operators.WeakDivergence()
    hgrad = Operators.Gradient()
    @. dY.soil.ϑ_l += -hdiv(-p.soil.K * hgrad(p.soil.ψ + model.coordinates.z))
    @. dY.soil.ρe_int +=
        -hdiv(
            -κ * hgrad(p.soil.T) -
            p.soil.K * ρe_int_l * hgrad(p.soil.ψ + model.coordinates.z),
        )
end

"""
    prognostic_vars(soil::SoilEnergyHydrology)

A function which returns the names of the prognostic variables
of `SoilEnergyHydrology`.
"""
ClimaLSM.prognostic_vars(soil::SoilEnergyHydrology) = (:ϑ_l, :θ_i, :ρe_int)

"""
    auxiliary_vars(soil::SoilEnergyHydrology)

A function which returns the names of the auxiliary variables
of `SoilEnergyHydrology`.
"""
ClimaLSM.auxiliary_vars(soil::SoilEnergyHydrology) = (:K, :ψ, :T)


"""
    make_update_aux(model::SoilEnergyHydrology)

An extension of the function `make_update_aux`, for the integrated
soil hydrology and energy model.

This function creates and returns a function which updates the auxiliary
variables `p.soil.variable` in place.

This has been written so as to work with Differential Equations.jl.
"""
function ClimaLSM.make_update_aux(model::SoilEnergyHydrology)
    function update_aux!(p, Y, t)
        @unpack ν, vg_α, vg_n, vg_m, Ksat, S_s, θ_r = model.hydrology_param_set
        @unpack ρc_s = model.thermal_param_set
        @. p.soil.K = hydraulic_conductivity(
            Ksat,
            vg_m,
            effective_saturation(ν, Y.soil.ϑ_l, θ_r),
        )
        @. p.soil.ψ = pressure_head(vg_α, vg_n, vg_m, θ_r, Y.soil.ϑ_l, ν, S_s)
        @. p.soil.T = temperature_from_ρe_int(Y.soil.ρe_int, Y.soil.θ_i, ρc_s, model.earth_param_set)
    end
    return update_aux!
end

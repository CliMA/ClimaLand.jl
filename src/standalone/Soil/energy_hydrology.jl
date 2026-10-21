"""
    EnergyHydrologyParameters{
            FT <: AbstractFloat,
            F <: Union{<:AbstractFloat, ClimaCore.Fields.Field},
            AP,
            C,
            PSE,
        }

A parameter structure for the integrated soil water and energy
 equation system.

- soil composition and retention model parameters defined
  in the interior
- an albedo parameterization
- Scalar parameters: emissivity, α, β, γ, γT_ref, Ω,
 roughness lengths z_0, d_ds ) (FT)

$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct EnergyHydrologyParameters{
    FT <: AbstractFloat,
    F <: Union{FT, ClimaCore.Fields.Field},
    AP,
    C,
    PSE,
}
    "The dry soil thermal conductivity, W/m/K"
    κ_dry::F
    "The saturated thermal conductivity of frozen soil, W/m/K"
    κ_sat_frozen::F
    "The saturated thermal conductivity of unfrozen soil, W/m/K"
    κ_sat_unfrozen::F
    "The volumetric heat capacity of dry soil, J/m^3/K (per volume dry soil, not per volume soil solids)"
    ρc_ds::F
    "The porosity of the soil (m^3/m^3)"
    ν::F
    "The volumetric fraction of the soil solids in organic matter (m^3/m^3)"
    ν_ss_om::F
    "The volumetric fraction of the soil solids in quartz (m^3/m^3)"
    ν_ss_quartz::F
    "The volumetric fraction of the soil solids in gravel (m^3/m^3)"
    ν_ss_gravel::F
    "The parameter α used in computing Kersten number, unitless"
    α::FT
    "The parameter β used in computing Kersten number, unitless"
    β::FT
    "The soil hydrology closure model: van Genuchten or Brooks and Corey"
    hydrology_cm::C
    "The saturated hydraulic conductivity (m/s)"
    K_sat::F
    "The specific storativity (1/m)"
    S_s::F
    "The residual water fraction (m^3/m^3)"
    θ_r::F
    "Ice impedance factor for the hydraulic conductivity"
    Ω::FT
    "Coefficient of viscosity factor for the hydraulic conductivity"
    γ::FT
    "Reference temperature for the viscosity factor"
    γT_ref::FT
    "Soil Albedo Parameterization"
    albedo::AP
    "Soil Emissivity"
    emissivity::FT
    "Roughness length for momentum"
    z_0m::FT
    "Roughness length for scalars"
    z_0b::FT
    "Maximum dry soil layer thickness under evaporation (m)"
    d_ds::FT
    "Exponent used in soil resistance to evaporation (unitless)"
    evap_p::FT
    "Multiplicative scalar used as fitting parameter in critical soil water content for evaporation (unitless)"
    evap_α::FT
    "Physical constants and clima-wide parameters"
    earth_param_set::PSE
end

Base.broadcastable(ps::EnergyHydrologyParameters) = tuple(ps)


"""
    EnergyHydrologyParameters(
        toml_dict;
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm,
        K_sat,
        S_s,
        θ_r,
        albedo = Soil.ConstantTwoBandSoilAlbedo{FT}(),
        kwargs...,)

EnergyHydrologyParameters has two constructors: float-type and toml dict based.
Additional parameters must be added manually: `ν`, `ν_ss_om`, `ν_ss_quartz`,
`ν_ss_gravel`, `hydrology_cm``, `K_sat`, `S_s`, and `θ_r`. All parameters can be
manually overridden via keyword arguments. Note, however, that certain parameters
must have the same type (e.g, if a field is supplied for porosity, it must be
supplied for all other parameters defined in the interior of the domain). Some
parameters are defined only on the surface of the domain (e.g albedo), while
other are defined everywhere (e.g. porosity). These are indicated with types `F`
and `SF`. If both dry/wet albedos and general albedos are given as keywords, the
dry/wet albedos will override the general albedos.

Please see the EnergyHydrologyParameters documentation for a complete list.
"""
function EnergyHydrologyParameters(
    toml_dict::CP.ParamDict;
    ν::F,
    ν_ss_om::F,
    ν_ss_quartz::F,
    ν_ss_gravel::F,
    hydrology_cm::C,
    K_sat::F,
    S_s::F,
    θ_r::F,
    albedo = Soil.ConstantTwoBandSoilAlbedo{CP.float_type(toml_dict)}(),
    emissivity = toml_dict["emissivity_bare_soil"],
    z_0m = toml_dict["soil_momentum_roughness_length"],
    z_0b = toml_dict["soil_scalar_roughness_length"],
    Ω = toml_dict["ice_impedance_omega"],
    d_ds = toml_dict["maximum_dry_soil_layer_depth"],
    evap_p = toml_dict["evaporation_exponent"],
    evap_α = toml_dict["evaporation_scalar"],
) where {F <: Union{<:AbstractFloat, ClimaCore.Fields.Field}, C}
    earth_param_set = LP.LandParameters(toml_dict)

    # Obtain parameters needed to calculate the derived parameters
    derived_param_name_map = (;
        :thermal_conductivity_of_quartz => :κ_quartz,
        :thermal_conductivity_of_soil_minerals => :κ_minerals,
        :thermal_conductivity_of_water_ice => :κ_ice,
        :thermal_conductivity_of_liquid_water => :κ_liq,
        :thermal_conductivity_of_organic_matter => :κ_om,
        :particle_density_quartz => :ρp_quartz,
        :particle_density_minerals => :ρp_minerals,
        :particle_density_organic_matter => :ρp_om,
        :vol_heat_capacity_quartz => :ρc_quartz,
        :vol_heat_capacity_organic_matter => :ρc_om,
        :vol_heat_capacity_minerals => :ρc_minerals,
    )
    p = CP.get_parameter_values(toml_dict, derived_param_name_map, "Land")
    # Particle density of the soil - per unit soil solids
    # Denoted ρ_ds in the Clima Design Docs (Equation 2.3)
    # where ν_ss_i = ν_i/(1-ν)
    ρp = @. (
        ν_ss_om * p.ρp_om +
        ν_ss_quartz * p.ρp_quartz +
        ν_ss_gravel * p.ρp_minerals +
        (1 - ν_ss_om - ν_ss_quartz - ν_ss_gravel) * p.ρp_minerals
    )
    # Volumetric heat capacity of soil solids - per unit volume soil solids
    # This is Equation 2.6a/2.10 in the Clima Design Docs
    # where ν_ss_i = ν_i/(1-ν)
    ρc_ss = @. (
        ν_ss_om * p.ρc_om +
        ν_ss_quartz * p.ρc_quartz +
        ν_ss_gravel * p.ρc_minerals +
        (1 - ν_ss_om - ν_ss_quartz - ν_ss_gravel) * p.ρc_minerals
    )
    # Volumetric heat capacity of dry soil - per unit volume of soil
    # This is denoted č_ds (Equation 2.11) and is used in equation 2.9
    ρc_ds = @. (1 - ν) * ρc_ss

    κ_solid =
        Soil.κ_solid.(ν_ss_om, ν_ss_quartz, p.κ_om, p.κ_quartz, p.κ_minerals)
    κ_dry = Soil.κ_dry.(ρp, ν, κ_solid, LP.K_therm(earth_param_set))
    κ_sat_frozen = Soil.κ_sat_frozen.(κ_solid, ν, p.κ_ice)
    κ_sat_unfrozen = Soil.κ_sat_unfrozen.(κ_solid, ν, p.κ_liq)

    name_map = (;
        :kersten_number_alpha => :α,
        :kersten_number_beta => :β,
        :temperature_factor_soil_hydraulic_conductivity => :γ,
        :temperature_reference_soil_hydraulic_conductivity => :γT_ref,
    )
    parameters = CP.get_parameter_values(toml_dict, name_map, "Land")
    PSE = typeof(earth_param_set)
    FT = CP.float_type(toml_dict)
    EnergyHydrologyParameters{FT, F, typeof(albedo), C, PSE}(;
        albedo,
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm,
        K_sat,
        S_s,
        θ_r,
        κ_dry,
        κ_sat_frozen,
        κ_sat_unfrozen,
        ρc_ds,
        earth_param_set,
        emissivity,
        z_0m,
        z_0b,
        Ω,
        d_ds,
        evap_p,
        evap_α,
        parameters...,
    )
end

"""
    EnergyHydrology <: AbstractSoilModel

A model for simulating the flow of water and heat
in a porous medium by solving the Richardson-Richards equation
and the heat equation, including terms for phase change.

A variety of boundary condition types are supported, including
`FluxBC`, `MoistureStateBC`/`TemperatureStateBC`,
`FreeDrainage` (only for the bottom of the domain),
and an `AtmosDrivenFluxBC` (under which radiative fluxes and
turbulent surface fluxes are computed and used as boundary conditions).
Please see the documentation for this boundary condition type for more
details.

# Fields
- `parameters`: The parameter sets
- `domain`: the soil domain, using ClimaCore.Domains
- `boundary_conditions`: the boundary conditions for RRE and heat, of type AbstractSoilBoundaryConditions
- `sources`: A tuple of sources, each of type AbstractSoilSource
- `lateral_flow`: A boolean flag which, when false, turns off the horizontal flow of water and heat
"""
mutable struct EnergyHydrology{FT, PS, D, BCRH, S} <: AbstractSoilModel{FT}
    const parameters::PS
    const domain::D
    const boundary_conditions::BCRH
    const sources::S
    const lateral_flow::Bool
end

"""
    EnergyHydrology{FT}(;
        parameters::PS
        domain::D,
        boundary_conditions::NamedTuple,
        sources::Tuple,
        lateral_flow::Bool = false,
    ) where {FT, D, PS}

A constructor for a `EnergyHydrology` model, which sets the default value
of the `lateral_flow` flag to false.
"""
function EnergyHydrology{FT}(;
    parameters::EnergyHydrologyParameters{FT, PSE},
    domain::D,
    boundary_conditions::NamedTuple,
    sources::Tuple,
    lateral_flow::Bool = false,
) where {FT, D, PSE}
    @assert !lateral_flow
    top_bc = boundary_conditions.top
    if typeof(top_bc) <: AtmosDrivenFluxBC
        # If the top BC indicates atmospheric conditions are driving the model
        # add baseflow as a sink term, add sublimation as a sink term
        subl_source =
            sublimation_source(Val(top_bc.prognostic_land_components), FT)
        subsurface_source = subsurface_runoff_source(top_bc.runoff)
        sources = append_source(subsurface_source, sources)
        sources = append_source(subl_source, sources)
    end
    args = (parameters, domain, boundary_conditions, sources)
    EnergyHydrology{FT, typeof.(args)...}(args..., lateral_flow)
end

function make_update_boundary_fluxes(model::EnergyHydrology)
    NVTX.@annotate function update_boundary_fluxes!(p, Y, t)
        Δz_top = model.domain.fields.Δz_top
        Δz_bottom = model.domain.fields.Δz_bottom
        soil_boundary_fluxes!(
            model.boundary_conditions.top,
            ClimaLand.TopBoundary(),
            model,
            Δz_top,
            Y,
            p,
            t,
        )

        soil_boundary_fluxes!(
            model.boundary_conditions.bottom,
            ClimaLand.BottomBoundary(),
            model,
            Δz_bottom,
            Y,
            p,
            t,
        )
    end
    return update_boundary_fluxes!
end



"""
    make_compute_exp_tendency(model::EnergyHydrology)

An extension of the function `make_compute_exp_tendency`, for the integrated
soil energy and heat equations, including phase change.

This function creates and returns a function which computes the entire
right hand side of the PDE for `Y.soil.ϑ_l, Y.soil.θ_i, Y.soil.ρe_int`,
and updates `dY.soil` in place with those values.
All of these quantities will be stepped explicitly.

This has been written so as to work with Differential Equations.jl.
"""
function ClimaLand.make_compute_exp_tendency(
    model::EnergyHydrology{FT},
) where {FT}
    NVTX.@annotate function compute_exp_tendency!(dY, Y, p, t)
        z = model.domain.fields.z

        dY.soil.∫F_vol_liq_water_dt .= 0
        dY.soil.∫F_e_dt .= 0

        # Don't update the prognostic variables we're stepping implicitly
        dY.soil.ϑ_l .= 0
        dY.soil.ρe_int .= 0

        # Note that soil ice content is only updated via source terms,
        # which are plus-added to dY.soil.θ_i, so we must zero it out here first.
        dY.soil.θ_i .= 0

        # Horizontal contributions
        horizontal_components!(
            dY,
            model.domain,
            Val(model.lateral_flow),
            model,
            p,
            z,
        )
        # Explicitly treated source terms
        # These change dY by +=, which is why we ".=" them above
        for src in model.sources
            if src.explicit
                ClimaLand.source!(dY, src, Y, p, model)
            end
        end
    end
    return compute_exp_tendency!
end

"""
    make_compute_imp_tendency(model::EnergyHydrology)

An extension of the function `make_compute_imp_tendency`, for the integrated
soil energy and heat equations, including phase change.

This version of this function computes the right hand side of the PDE for
`Y.soil.ϑ_l`, which is the only quantity we currently step implicitly.

This has been written so as to work with Differential Equations.jl.
"""
function ClimaLand.make_compute_imp_tendency(
    model::EnergyHydrology{FT},
) where {FT}
    NVTX.@annotate function compute_imp_tendency!(dY, Y, p, t)
        z = model.domain.fields.z
        rre_top_flux_bc = p.soil.top_bc.water
        rre_bottom_flux_bc = p.soil.bottom_bc.water
        @. dY.soil.∫F_vol_liq_water_dt = -(rre_top_flux_bc - rre_bottom_flux_bc) # These fluxes appear in implicit terms, we step them implicitly
        heat_top_flux_bc = p.soil.top_bc.heat
        heat_bottom_flux_bc = p.soil.bottom_bc.heat
        @. dY.soil.∫F_e_dt = -(heat_top_flux_bc - heat_bottom_flux_bc) # These fluxes appear in implicit terms, we step them implicitly
        interpc2f = Operators.InterpolateC2F()
        gradc2f = Operators.GradientC2F()

        # Without topography only
        # In Cartesian coordinates, W (z^) = Cov3 (z^)= Contra3 (n^ = z^)
        # In spherical coordinates, W (r^) = Cov3 (r^) = Contra3 (n^ = r^)
        # Passing WVector to gradient BC is passing a normal flux.

        # Richards-Richardson RHS
        @. p.soil.top_bc_wvec = Geometry.WVector(rre_top_flux_bc)
        @. p.soil.bottom_bc_wvec = Geometry.WVector(rre_bottom_flux_bc)
        divf2c_rre = Operators.DivergenceF2C(
            top = Operators.SetValue(p.soil.top_bc_wvec),
            bottom = Operators.SetValue(p.soil.bottom_bc_wvec),
        )
        # GradC2F returns a Covariant3Vector, so no need to convert.
        @. dY.soil.ϑ_l =
            -(divf2c_rre(-interpc2f(p.soil.K) * gradc2f(p.soil.ψ + z)))

        # Heat equation RHS
        # Reuse the same scratch space:
        @. p.soil.top_bc_wvec = Geometry.WVector(heat_top_flux_bc)
        @. p.soil.bottom_bc_wvec = Geometry.WVector(heat_bottom_flux_bc)
        divf2c_heat = Operators.DivergenceF2C(
            top = Operators.SetValue(p.soil.top_bc_wvec),
            bottom = Operators.SetValue(p.soil.bottom_bc_wvec),
        )

        # GradC2F returns a Covariant3Vector, so no need to convert.
        @. dY.soil.ρe_int =
            -divf2c_heat(
                -interpc2f(p.soil.κ) * gradc2f(p.soil.T) -
                interpc2f(
                    volumetric_internal_energy_liq(
                        p.soil.T,
                        model.parameters.earth_param_set,
                    ) * p.soil.K,
                ) * gradc2f(p.soil.ψ + z),
            )

        @. dY.soil.θ_i = 0

        # Source terms
        # These change dY by +=, which is why we ".=" above
        for src in model.sources
            if !src.explicit
                ClimaLand.source!(dY, src, Y, p, model)
            end
        end
    end
    return compute_imp_tendency!
end

function ClimaLand.make_update_implicit_aux(model::EnergyHydrology)
    NVTX.@annotate function update_imp_aux!(p, Y, t)
        (; ν, hydrology_cm, S_s, θ_r, ρc_ds, earth_param_set) = model.parameters
        @. p.soil.θ_l =
            volumetric_liquid_fraction(Y.soil.ϑ_l, ν - Y.soil.θ_i, θ_r)
        @. p.soil.T = temperature_from_ρe_int(
            Y.soil.ρe_int,
            Y.soil.θ_i,
            volumetric_heat_capacity(
                p.soil.θ_l,
                Y.soil.θ_i,
                ρc_ds,
                earth_param_set,
            ),
            earth_param_set,
        )
        @. p.soil.ψ =
            pressure_head(hydrology_cm, θ_r, Y.soil.ϑ_l, ν - Y.soil.θ_i, S_s)
    end
    return update_imp_aux!
end

function ClimaLand.make_update_implicit_boundary_fluxes(model::EnergyHydrology)
    ubf! = make_update_boundary_fluxes(model)
    NVTX.@annotate function update_imp_bf!(p, Y, t)
        if haskey(p.soil, :dfluxBCdY)
            ubf!(p, Y, t)
        end
    end
    return update_imp_bf!
end

"""
    ClimaLand.make_compute_jacobian(model::EnergyHydrology{FT}) where {FT}

Creates and returns the compute_jacobian! function for the EnergyHydrology model.
This updates the contribution for the soil liquid water content only.

Using this Jacobian with a backwards Euler timestepper is equivalent
to using the modified Picard scheme of Celia et al. (1990).
"""
function ClimaLand.make_compute_jacobian(model::EnergyHydrology{FT}) where {FT}
    NVTX.@annotate function compute_jacobian!(
        jacobian::MatrixFields.FieldMatrixWithSolver,
        Y,
        p,
        dtγ,
        t,
    )
        (; matrix) = jacobian
        (; ν, hydrology_cm, S_s, θ_r, ρc_ds, earth_param_set) = model.parameters

        # Create divergence operator
        divf2c_op = Operators.DivergenceF2C()
        divf2c_matrix = MatrixFields.operator_matrix(divf2c_op)
        # Create gradient operator, and set gradient at boundaries to 0
        gradc2f_op = Operators.GradientC2F(
            top = Operators.SetGradient(Geometry.WVector(FT(0))),
            bottom = Operators.SetGradient(Geometry.WVector(FT(0))),
        )
        gradc2f_matrix = MatrixFields.operator_matrix(gradc2f_op)
        # Create interpolation operator
        interpc2f_op = Operators.InterpolateC2F(
            bottom = Operators.Extrapolate(),
            top = Operators.Extrapolate(),
        )

        # The derivative of the residual with respect to the prognostic variable
        ∂ϑres∂ϑ = matrix[@name(soil.ϑ_l), @name(soil.ϑ_l)]
        ∂ρeres∂ρe = matrix[@name(soil.ρe_int), @name(soil.ρe_int)]
        ∂ρeres∂ϑ = matrix[@name(soil.ρe_int), @name(soil.ϑ_l)]

        # Derivatives with respect to ϑ:

        # Precompute intermediate quantities to improve performance
        # due to fusing of broadcasted expressions involving matrices
        # First, the gradient of ∂ψ∂ϑ
        # This term is used again below, so we do not alter it once we have made it
        @. p.soil.bidiag_matrix_scratch =
            gradc2f_matrix() * MatrixFields.DiagonalMatrixRow(
                ClimaLand.Soil.dψdϑ(
                    hydrology_cm,
                    Y.soil.ϑ_l,
                    ν - Y.soil.θ_i, #ν_eff
                    θ_r,
                    S_s,
                ),
            )
        # Now the full Darcy flux term. This term is the one that gets altered with the BC
        # contribution in place, below
        @. p.soil.full_bidiag_matrix_scratch =
            MatrixFields.DiagonalMatrixRow(interpc2f_op(-p.soil.K)) *
            p.soil.bidiag_matrix_scratch

        # If the top BC is a `MoistureStateBC`, add the term from the top BC
        #  flux before applying divergence
        if haskey(p.soil, :dfluxBCdY)
            dfluxBCdY = p.soil.dfluxBCdY
            topBC_op = Operators.SetBoundaryOperator(
                top = Operators.SetValue(dfluxBCdY),
                bottom = Operators.SetValue(
                    Geometry.Covariant3Vector(zero(FT)),
                ),
            )
            # Note: need to pass 3D field on faces to `topBC_op`. Interpolating `K` to faces
            #  for this is inefficient - we should find a better solution.
            @. p.soil.topBC_scratch = topBC_op(
                Geometry.Covariant3Vector(zero(interpc2f_op(p.soil.K))),
            )
            @. p.soil.full_bidiag_matrix_scratch +=
                MatrixFields.LowerDiagonalMatrixRow(p.soil.topBC_scratch)
        end
        # dtγ can be an ITime or a float
        negative_dtγ = FT(-float(dtγ))
        @. ∂ϑres∂ϑ =
            negative_dtγ *
            (divf2c_matrix() * p.soil.full_bidiag_matrix_scratch) - (I,)

        # Now create the flux term for ∂ρe∂ϑ using bidiag_matrix_scratch
        # This overwrites full_bidiag_matrix_scratch
        @. p.soil.full_bidiag_matrix_scratch =
            MatrixFields.DiagonalMatrixRow(
                -interpc2f_op(
                    volumetric_internal_energy_liq(
                        p.soil.T,
                        model.parameters.earth_param_set,
                    ) * p.soil.K,
                ),
            ) * p.soil.bidiag_matrix_scratch
        @. ∂ρeres∂ϑ =
            negative_dtγ * (divf2c_matrix() * p.soil.full_bidiag_matrix_scratch)

        # Now overwrite bidiag_matrix_scratch and full_bidiag scratch for the ρe ρe bidiagonal
        @. p.soil.bidiag_matrix_scratch =
            gradc2f_matrix() * MatrixFields.DiagonalMatrixRow(
                1 / ClimaLand.Soil.volumetric_heat_capacity(
                    p.soil.θ_l,
                    Y.soil.θ_i,
                    ρc_ds,
                    earth_param_set,
                ),
            )
        @. p.soil.full_bidiag_matrix_scratch =
            MatrixFields.DiagonalMatrixRow(interpc2f_op(-p.soil.κ)) *
            p.soil.bidiag_matrix_scratch
        @. ∂ρeres∂ρe =
            negative_dtγ *
            (divf2c_matrix() * p.soil.full_bidiag_matrix_scratch) - (I,)
    end
    return compute_jacobian!
end

"""
   horizontal_components!(dY::ClimaCore.Fields.FieldVector,
                          domain::Union{HybridBox, SphericalShell},
                          lateral_flow::Val{true},
                          model::EnergyHydrology,
                          p::NamedTuple)

Updates dY in place by adding in the tendency terms resulting from
horizontal derivative operators for the `EnergyHydrology` model,
 in the case of a hybrid box or
spherical shell domain with the model
`lateral_flag` set to true.

The horizontal contributions are
computed using the WeakDivergence and Gradient operators.
"""
NVTX.@annotate function horizontal_components!(
    dY::ClimaCore.Fields.FieldVector,
    domain::Union{HybridBox, SphericalShell},
    lateral_flow::Val{true},
    model::EnergyHydrology,
    p::NamedTuple,
    z::ClimaCore.Fields.Field,
)
    hdiv = Operators.Divergence{Operators.WeakForm}()
    hgrad = Operators.Gradient()
    # The flux is already covariant, from hgrad, so no need to convert.
    @. dY.soil.ϑ_l += -hdiv(-p.soil.K * hgrad(p.soil.ψ + z))
    @. dY.soil.ρe_int +=
        -hdiv(
            -p.soil.κ * hgrad(p.soil.T) -
            p.soil.K *
            volumetric_internal_energy_liq(
                p.soil.T,
                model.parameters.earth_param_set,
            ) *
            hgrad(p.soil.ψ + z),
        )
end

"""
    prognostic_vars(soil::EnergyHydrology)

A function which returns the names of the prognostic variables
of `EnergyHydrology`.
"""
ClimaLand.prognostic_vars(soil::EnergyHydrology) =
    (:ϑ_l, :θ_i, :ρe_int, :∫F_vol_liq_water_dt, :∫F_e_dt)

"""
    prognostic_types(soil::EnergyHydrology{FT}) where {FT}

A function which returns the types of the prognostic variables
of `EnergyHydrology`.
"""
ClimaLand.prognostic_types(soil::EnergyHydrology{FT}) where {FT} =
    (FT, FT, FT, FT, FT)

ClimaLand.prognostic_domain_names(soil::EnergyHydrology) =
    (:subsurface, :subsurface, :subsurface, :surface, :surface)
"""
    auxiliary_vars(soil::EnergyHydrology)

A function which returns the names of the auxiliary variables
of `EnergyHydrology`.
"""
ClimaLand.auxiliary_vars(soil::EnergyHydrology) = (
    :total_water,
    :total_energy,
    :K,
    :ψ,
    :θ_l,
    :T,
    :κ,
    :Tf_depressed,
    :bidiag_matrix_scratch,
    :full_bidiag_matrix_scratch,
    boundary_vars(soil.boundary_conditions.top, ClimaLand.TopBoundary())...,
    boundary_vars(
        soil.boundary_conditions.bottom,
        ClimaLand.BottomBoundary(),
    )...,
)

"""
    auxiliary_types(soil::EnergyHydrology{FT}) where {FT}

A function which returns the types of the auxiliary variables
of `EnergyHydrology`.
"""
ClimaLand.auxiliary_types(soil::EnergyHydrology{FT}) where {FT} = (
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    FT,
    MatrixFields.BidiagonalMatrixRow{Geometry.Covariant3Vector{FT}},
    MatrixFields.BidiagonalMatrixRow{Geometry.Covariant3Vector{FT}},
    boundary_var_types(
        soil,
        soil.boundary_conditions.top,
        ClimaLand.TopBoundary(),
    )...,
    boundary_var_types(
        soil,
        soil.boundary_conditions.bottom,
        ClimaLand.BottomBoundary(),
    )...,
)

ClimaLand.auxiliary_domain_names(soil::EnergyHydrology) = (
    :surface,
    :surface,
    :subsurface,
    :subsurface,
    :subsurface,
    :subsurface,
    :subsurface,
    :subsurface,
    :subsurface_face,
    :subsurface_face,
    boundary_var_domain_names(
        soil.boundary_conditions.top,
        ClimaLand.TopBoundary(),
    )...,
    boundary_var_domain_names(
        soil.boundary_conditions.bottom,
        ClimaLand.BottomBoundary(),
    )...,
)
"""
    make_update_aux(model::EnergyHydrology)

An extension of the function `make_update_aux`, for the integrated
soil hydrology and energy model.

This function creates and returns a function which updates the auxiliary
variables `p.soil.variable` in place.

This has been written so as to work with Differential Equations.jl.
"""
function ClimaLand.make_update_aux(model::EnergyHydrology)
    NVTX.@annotate function update_aux!(p, Y, t)
        (;
            ν,
            hydrology_cm,
            K_sat,
            S_s,
            θ_r,
            Ω,
            γ,
            α,
            β,
            ν_ss_om,
            ν_ss_gravel,
            ν_ss_quartz,
            γT_ref,
            κ_sat_frozen,
            κ_sat_unfrozen,
            ρc_ds,
            earth_param_set,
            albedo,
        ) = model.parameters

        @. p.soil.θ_l =
            volumetric_liquid_fraction(Y.soil.ϑ_l, ν - Y.soil.θ_i, θ_r)

        update_albedo!(
            model.boundary_conditions.top,
            albedo,
            p,
            model.domain,
            model.parameters,
        )

        @. p.soil.κ = thermal_conductivity(
            model.parameters.κ_dry,
            kersten_number(
                Y.soil.θ_i,
                relative_saturation(p.soil.θ_l, Y.soil.θ_i, ν),
                α,
                β,
                ν_ss_om,
                ν_ss_quartz,
                ν_ss_gravel,
            ),
            κ_sat(p.soil.θ_l, Y.soil.θ_i, κ_sat_unfrozen, κ_sat_frozen),
        )

        @. p.soil.T = temperature_from_ρe_int(
            Y.soil.ρe_int,
            Y.soil.θ_i,
            volumetric_heat_capacity(
                p.soil.θ_l,
                Y.soil.θ_i,
                ρc_ds,
                earth_param_set,
            ),
            earth_param_set,
        )

        @. p.soil.K =
            impedance_factor(Y.soil.θ_i / (p.soil.θ_l + Y.soil.θ_i - θ_r), Ω) *
            viscosity_factor(p.soil.T, γ, γT_ref) *
            hydraulic_conductivity(
                hydrology_cm,
                K_sat,
                effective_saturation(ν, Y.soil.ϑ_l, θ_r),
            )
        @. p.soil.ψ =
            pressure_head(hydrology_cm, θ_r, Y.soil.ϑ_l, ν - Y.soil.θ_i, S_s)
        _ρ_l = LP.ρ_cloud_liq(earth_param_set)
        _ρ_i = LP.ρ_cloud_ice(earth_param_set)
        _LH_f0 = LP.LH_f0(earth_param_set)
        _T_freeze = LP.T_freeze(earth_param_set)
        _grav = LP.grav(earth_param_set)
        @. p.soil.Tf_depressed = soil_Tf_depressed(
            p.soil.θ_l,
            Y.soil.θ_i,
            ν,
            θ_r,
            hydrology_cm,
            _ρ_l,
            _ρ_i,
            _T_freeze,
            _grav,
            _LH_f0,
        )

        total_liq_water_vol_per_area!(p.soil.total_water, model, Y, p, t)
        total_energy_per_area!(p.soil.total_energy, model, Y, p, t)
    end
    return update_aux!
end

"""
    PhaseChange{FT} <: AbstractSoilSource{FT}

PhaseChange source type; treated explicitly in all
prognostic variables.
"""
@kwdef struct PhaseChange{FT} <: AbstractSoilSource{FT}
    explicit::Bool = true
end

"""
    source!(dY::ClimaCore.Fields.FieldVector,
            src::PhaseChange{FT},
            Y::ClimaCore.Fields.FieldVector,
            p::NamedTuple,
            model
            )

Computes the source terms for phase change
explicitly in time.

"""
NVTX.@annotate function ClimaLand.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::PhaseChange,
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    model,
)
    params = model.parameters
    (; ν, ρc_ds, θ_r, hydrology_cm, earth_param_set) = params
    _ρ_l = LP.ρ_cloud_liq(earth_param_set)
    _ρ_i = LP.ρ_cloud_ice(earth_param_set)
    Δz = model.domain.fields.Δz # center face distance
    _LH_f0 = LP.LH_f0(earth_param_set)
    _T_freeze = LP.T_freeze(earth_param_set)
    _grav = LP.grav(earth_param_set)

    # Pass in physical constants as floats, to avoid issue broadcasting over fields and structs
    @. dY.soil.ϑ_l +=
        -phase_change_source(
            p.soil.θ_l,
            Y.soil.θ_i,
            p.soil.T,
            thermal_time(
                volumetric_heat_capacity(
                    p.soil.θ_l,
                    Y.soil.θ_i,
                    ρc_ds,
                    earth_param_set,
                ),
                Δz,
                p.soil.κ,
            ),
            ν,
            θ_r,
            hydrology_cm,
            _ρ_i,
            _ρ_l,
            _LH_f0,
            _T_freeze,
            _grav,
        )
    @. dY.soil.θ_i +=
        (_ρ_l / _ρ_i) * phase_change_source(
            p.soil.θ_l,
            Y.soil.θ_i,
            p.soil.T,
            thermal_time(
                volumetric_heat_capacity(
                    p.soil.θ_l,
                    Y.soil.θ_i,
                    ρc_ds,
                    earth_param_set,
                ),
                Δz,
                p.soil.κ,
            ),
            ν,
            θ_r,
            hydrology_cm,
            _ρ_i,
            _ρ_l,
            _LH_f0,
            _T_freeze,
            _grav,
        )

    # These source/sink terms conserve total water, so we do not include them when checking conservation because their
    # contribution is zero by design.
end


"""
    SoilSublimation{FT} <: AbstractSoilSource{FT}

Soil Sublimation source type. Treated explicitly since the turbulent
fluxes for SHF and LHF are also treated explicitly.
"""
@kwdef struct SoilSublimation{FT} <: AbstractSoilSource{FT}
    explicit::Bool = true
end
"""
    source!(dY::ClimaCore.Fields.FieldVector,
            src::SoilSublimation{FT},
            Y::ClimaCore.Fields.FieldVector,
            p::NamedTuple,
            model
            )

Updates dY.soil.θ_i in place with a term due to sublimation; this only affects
the surface layer of soil.

"""
NVTX.@annotate function ClimaLand.source!(
    dY::ClimaCore.Fields.FieldVector,
    src::SoilSublimation{FT},
    Y::ClimaCore.Fields.FieldVector,
    p::NamedTuple,
    model,
) where {FT}
    _ρ_i = FT(LP.ρ_cloud_ice(model.parameters.earth_param_set))
    _ρ_l = FT(LP.ρ_cloud_liq(model.parameters.earth_param_set))
    z = model.domain.fields.z
    Δz_top = model.domain.fields.Δz_top # this returns the center-face distance, not layer thickness
    @. dY.soil.θ_i +=
        -p.soil.turbulent_fluxes.vapor_flux_ice * _ρ_l / _ρ_i *
        heaviside(z + 2 * Δz_top) / (2 * Δz_top) # only apply to top layer, recall that z is negative
    @. dY.soil.∫F_vol_liq_water_dt += -p.soil.turbulent_fluxes.vapor_flux_ice # The integral of the source is designed to be this
end

## The functions below are required to be defined
## for use with the `AtmosDrivenFluxBC` upper
## boundary conditions. They return the
## surface conditions necessary for computing
## radiative, sensible and latent heat fluxes
## as well as evaporation.

"""
    ClimaLand.component_temperature(
        model::EnergyHydrology{FT},
        Y,
        p,
    ) where {FT}

Returns the surface temperature field of the
`EnergyHydrology` soil model.

For atmospherically driven soil (`AtmosDrivenFluxBC`), this is the soil skin
temperature `p.soil.turbulent_fluxes.T_sfc`, solved for from the surface energy
balance (see `update_soil_surface_temperature!`); otherwise it is the
temperature at the center of the first soil layer.
"""
function ClimaLand.component_temperature(model::EnergyHydrology, Y, p)
    return soil_surface_temperature(model.boundary_conditions.top, p)
end

"""
    ClimaLand.surface_emissivity(
        model::EnergyHydrology,
        Y,
        p,
    )

Returns the surface emissivity field of the
`EnergyHydrology` soil model.
"""
function ClimaLand.surface_emissivity(model::EnergyHydrology, Y, p)
    return model.parameters.emissivity
end

"""
    ClimaLand.component_specific_humidity(
        model::EnergyHydrology,
        Y,
        p,
        forcing = (;
            h_atmos = model.boundary_conditions.top.atmos.h,
            T_atmos = p.drivers.T,
            q_atmos = p.drivers.q,
            reference = model.boundary_conditions.top.atmos.reference
        ),
    )

Return the surface specific humidity for the `EnergyHydrology` soil model,
updating `p.soil.q_sfc` in place. The surface air density is extrapolated from
the air temperature `T_atmos` and specific humidity `q_atmos` at the reference
height `h_atmos` of `forcing`, which default to the atmospheric forcing; below
a canopy, the sub-canopy forcing is passed. Please be aware that the reference level
may be different between the two cases.
"""
function ClimaLand.component_specific_humidity(
    model::EnergyHydrology,
    Y,
    p,
    forcing = (;
        h_atmos = model.boundary_conditions.top.atmos.h,
        T_atmos = p.drivers.T,
        q_atmos = p.drivers.q,
        reference = model.boundary_conditions.top.atmos.reference,
    ),
)
    (; h_atmos, T_atmos, q_atmos, reference) = forcing
    earth_param_set = get_earth_param_set(model)
    surface_flux_params =
        LP.surface_fluxes_parameters(model.parameters.earth_param_set)
    T_sfc = component_temperature(model, Y, p)
    Tf_depressed_sfc =
        ClimaLand.Domains.top_center_to_surface(p.soil.Tf_depressed)
    ψ_sfc = ClimaLand.Domains.top_center_to_surface(p.soil.ψ)
    z_sfc = ClimaLand.surface_elevation(model, Y, p)
    displ = ClimaLand.surface_displacement_height(model, Y, p)
    z_0m = model.parameters.z_0m
    ρ_sfc = @.lazy(
        ClimaLand.compute_ρ_sfc(
            surface_flux_params,
            T_atmos,
            p.drivers.P,
            q_atmos,
            ClimaLand.compute_Δz_eff(h_atmos, z_0m, displ, reference),
            T_sfc,
        ),
    )
    @. p.soil.q_sfc = soil_specific_humidity(
        T_sfc,
        ρ_sfc,
        ψ_sfc,
        Tf_depressed_sfc,
        earth_param_set,
    )
    return p.soil.q_sfc
end

function soil_specific_humidity(
    T_sfc,
    ρ_sfc,
    ψ_sfc,
    Tf_depressed_sfc,
    earth_param_set,
)
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    g = LP.grav(earth_param_set)
    R = LP.gas_constant(earth_param_set)
    M_w = LP.molar_mass_water(earth_param_set)
    # Compute q_soil using ice or liquid as appropriate
    if T_sfc > Tf_depressed_sfc # liquid water evaporation
        q_sfc =
            Thermodynamics.q_vap_saturation(
                thermo_params,
                T_sfc,
                ρ_sfc,
                Thermodynamics.Liquid(),
            ) * exp(g * min(ψ_sfc, zero(ψ_sfc)) * M_w / (R * T_sfc))
    else
        q_sfc = Thermodynamics.q_vap_saturation(
            thermo_params,
            T_sfc,
            ρ_sfc,
            Thermodynamics.Ice(),
        )
    end
    return q_sfc
end

function ClimaLand.get_update_surface_humidity_function(
    model::EnergyHydrology,
    Y,
    p,
)
    function update_q_vap_sfc_at_a_point(
        ζ,
        param_set,
        thermo_params,
        inputs,
        scheme,
        T_sfc,
        u_star,
        z_0m,
        z_0b,
        g_liq,
        β_ice,
        Tf_depressed,
        qsat_sfc,
    )
        g_h = SurfaceFluxes.heat_conductance(
            param_set,
            ζ,
            u_star,
            inputs,
            z_0m,
            z_0b,
            scheme,
        )
        q_air = SurfaceFluxes.interior_vapor_specific_humidity(inputs)
        # Sublimation; at T_sfc = Tf_depressed, qsat_sfc is over ice
        frozen = inputs.T_sfc_guess <= Tf_depressed
        w = soil_evaporation_beta(q_air, qsat_sfc, g_liq, g_h, β_ice, frozen)
        return w * qsat_sfc + (1 - w) * q_air
    end
    # Closure; p.soil.q_sfc holds the saturated humidity from turbulent_fluxes!
    qsat_sfc = p.soil.q_sfc
    Tf_depressed_sfc =
        ClimaLand.Domains.top_center_to_surface(p.soil.Tf_depressed)
    ν_sfc = ClimaLand.Domains.top_center_to_surface(model.parameters.ν)
    θ_i_sfc = ClimaLand.Domains.top_center_to_surface(Y.soil.θ_i)
    g_soil_sfc =
        soil_surface_vapor_conductance!(p.soil.sfc_scratch, model, Y, p)
    update_q_vap_sfc_field(g_liq, β_ice, Tf_depressed, qsat_sfc) =
        (args...) -> update_q_vap_sfc_at_a_point(
            args...,
            g_liq,
            β_ice,
            Tf_depressed,
            qsat_sfc,
        )
    return @. lazy(
        update_q_vap_sfc_field(
            g_soil_sfc,
            frozen_soil_vapor_weight(θ_i_sfc, ν_sfc),
            Tf_depressed_sfc,
            qsat_sfc,
        ),
    )
end

"""
    ClimaLand.get_∂T_sfc∂T_function(model::EnergyHydrology, Y, p)
    ClimaLand.get_∂q_sfc∂T_function(model::EnergyHydrology, Y, p)

Return `nothing`: `p.soil.turbulent_fluxes` stores no derivatives of the
fluxes with respect to the soil temperature, so
`ClimaLand.turbulent_fluxes_from_output` omits them from the skin temperature
solve.
"""
ClimaLand.get_∂T_sfc∂T_function(model::EnergyHydrology, Y, p) = nothing
ClimaLand.get_∂q_sfc∂T_function(model::EnergyHydrology, Y, p) = nothing

function ClimaLand.surface_roughness_model(
    model::EnergyHydrology{FT},
    Y,
    p,
) where {FT}
    params = model.parameters
    return SurfaceFluxes.ConstantRoughnessParams{FT}(params.z_0m, params.z_0b)
end


function ClimaLand.get_drivers(model::EnergyHydrology)
    bc = model.boundary_conditions.top
    if typeof(bc) <: AtmosDrivenFluxBC{
        <:PrescribedAtmosphere,
        <:AbstractRadiativeDrivers,
        <:AbstractRunoffModel,
    }
        return (bc.atmos, bc.radiation)
    else
        return ()
    end
end

"""
    ClimaLand.total_liq_water_vol_per_area!(
        surface_field,
        model::EnergyHydrology,
        Y,
        p,
        t,
    )

A function which updates `surface_field` in place with the value for
the total liquid water volume per unit ground area for the `EnergyHydrology`.

The water in all phases is accounted for by converting ice volume to liquid water
volume using the ratio of the density of ice to the density of water.
"""
function ClimaLand.total_liq_water_vol_per_area!(
    surface_field,
    model::EnergyHydrology,
    Y,
    p,
    t,
)
    earth_param_set = model.parameters.earth_param_set
    _ρ_liq = LP.ρ_cloud_liq(earth_param_set)
    _ρ_ice = LP.ρ_cloud_ice(earth_param_set)
    ClimaCore.Operators.column_integral_definite!(
        surface_field,
        Y.soil.ϑ_l .+ Y.soil.θ_i .* _ρ_ice ./ _ρ_liq, # allocates
    )
    return nothing
end

"""
    ClimaLand.total_energy_per_area!(
        surface_field,
        model::EnergyHydrology,
        Y,
        p,
        t,
    )

A function which updates `surface_field` in place with the value for
the total energy per unit ground area for the `EnergyHydrology`.
"""
function ClimaLand.total_energy_per_area!(
    surface_field,
    model::EnergyHydrology,
    Y,
    p,
    t,
)
    ClimaCore.Operators.column_integral_definite!(surface_field, Y.soil.ρe_int)
    return nothing
end

"""
    soil_tortuosity(ν, θ_r, θ_i)

Compute the tortuosity factor `θ_a^2.5 / ν` [-] for water vapor diffusion
through the dry surface layer of the soil, following Equation (1) of Shokri,
Lehmann, and Or (2008), Geophys. Res. Lett., 35, L19407,
doi:10.1029/2008GL035230, with porosity `ν`, residual water fraction `θ_r`,
and ice fraction `θ_i`.

The vapor diffuses through the air-filled pore space of the dry layer, where
the liquid water content is residual, so `θ_a = ν - θ_r - θ_i`, independent of
the moisture of the soil below the dry layer (as in Swenson and Lawrence 2014,
J. Geophys. Res. Atmos., 119, 10299–10312, doi:10.1002/2014JD022314, and the
Community Land Model).
"""
function soil_tortuosity(ν, θ_r, θ_i)
    θ_a = ν - θ_r - θ_i
    safe_θ_a = max(θ_a, eps(typeof(θ_a)))
    return safe_θ_a^oftype(safe_θ_a, 2.5) / ν
end

"""
    soil_conductance(S_l, S_c, d_ds, p, α, _D_vapor, ν, θ_r, θ_i)

Computes the conductance of the top of the soil column to
water vapor diffusion, as a function of the surface
effective liquid water saturation `S_l`, critical saturation `S_c`,
other soil parameters, and diffusivity of vapor in air.

The conductance is the inverse of the diffusive resistance
`dsl / (D_vapor τ_a)` of the dry soil layer of thickness `dsl`, with the
tortuosity factor `τ_a` of the dry layer (see [`soil_tortuosity`](@ref)).
When no dry layer has formed (`dsl = 0`), the resistance vanishes and the
conductance is unbounded.
"""
function soil_conductance(S_l, S_c, d_ds, p, α, _D_vapor, ν, θ_r, θ_i)
    dsl = dry_soil_layer_thickness(S_l, α * S_c, d_ds, p)
    τ_a = soil_tortuosity(ν, θ_r, θ_i)
    r_soil = dsl / (_D_vapor * τ_a) # [s/m]
    return inv(max(r_soil, eps(typeof(r_soil)))) # [m/s]
end

"""
    dry_soil_layer_thickness(S_l, αS_c, d_ds, p)

Returns the maximum dry soil layer thickness that can develop under vapor flux; 
this is used when computing the soil resistance to vapor flux similar to
Swenson et al (2012)/Sakaguchi and Zeng (2009).

Note that S_l ∈ (0,1].
"""
function dry_soil_layer_thickness(S_l, αS_c, d_ds, p)
    αS_c = clamp(αS_c, zero(αS_c), one(αS_c))
    return S_l < αS_c ? d_ds * ((αS_c - S_l) / S_l)^p :
           zero(d_ds * (αS_c - S_l) / S_l * p)
end

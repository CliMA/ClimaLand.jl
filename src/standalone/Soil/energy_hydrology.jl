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

$(DocStringExtensions.FIELDS)
"""
struct EnergyHydrology{FT, PS, D, BCRH, S} <: AbstractSoilModel{FT}
    "The parameter sets"
    parameters::PS
    "the soil domain, using ClimaCore.Domains"
    domain::D
    "the boundary conditions for RRE and heat, of type AbstractSoilBoundaryConditions"
    boundary_conditions::BCRH
    "A tuple of sources, each of type AbstractSoilSource"
    sources::S
    "A boolean flag which, when false, turns off the horizontal flow of water and heat"
    lateral_flow::Bool
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
    function update_boundary_fluxes!(p, Y, t)
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
    function compute_exp_tendency!(dY, Y, p, t)
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
    function compute_imp_tendency!(dY, Y, p, t)
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
    function update_imp_aux!(p, Y, t)
        (; ν, hydrology_cm, S_s, θ_r, ρc_ds, earth_param_set) = model.parameters
        @. p.soil.T = temperature_from_ρe_int(
            Y.soil.ρe_int,
            Y.soil.θ_i,
            volumetric_heat_capacity(
                min(ν - Y.soil.θ_i, Y.soil.ϑ_l), # compute θ_l
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
    function update_imp_bf!(p, Y, t)
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
    function compute_jacobian!(
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
            gradc2f_matrix() ⋅ MatrixFields.DiagonalMatrixRow(
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
            MatrixFields.DiagonalMatrixRow(interpc2f_op(-p.soil.K)) ⋅
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
        @. ∂ϑres∂ϑ =
            FT(-1) *
            float(dtγ) *
            (divf2c_matrix() ⋅ p.soil.full_bidiag_matrix_scratch) - (I,)

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
            ) ⋅ p.soil.bidiag_matrix_scratch
        @. ∂ρeres∂ϑ =
            FT(-1) *
            float(dtγ) *
            (divf2c_matrix() ⋅ p.soil.full_bidiag_matrix_scratch) - (I,)

        # Now overwrite bidiag_matrix_scratch and full_bidiag scratch for the ρe ρe bidiagonal
        @. p.soil.bidiag_matrix_scratch =
            gradc2f_matrix() ⋅ MatrixFields.DiagonalMatrixRow(
                1 / ClimaLand.Soil.volumetric_heat_capacity(
                    p.soil.θ_l,
                    Y.soil.θ_i,
                    ρc_ds,
                    earth_param_set,
                ),
            )
        @. p.soil.full_bidiag_matrix_scratch =
            MatrixFields.DiagonalMatrixRow(interpc2f_op(-p.soil.κ)) ⋅
            p.soil.bidiag_matrix_scratch
        @. ∂ρeres∂ρe =
            FT(-1) *
            float(dtγ) *
            (divf2c_matrix() ⋅ p.soil.full_bidiag_matrix_scratch) - (I,)
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
function horizontal_components!(
    dY::ClimaCore.Fields.FieldVector,
    domain::Union{HybridBox, SphericalShell},
    lateral_flow::Val{true},
    model::EnergyHydrology,
    p::NamedTuple,
    z::ClimaCore.Fields.Field,
)
    hdiv = Operators.WeakDivergence()
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
    function update_aux!(p, Y, t)
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
function ClimaLand.source!(
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

Soil Sublimation source type. Used to defined a method
of `ClimaLand.source!` for soil sublimation; treated implicitly
in ϑ_l, ρe_int but explicitly in θ_i.
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
function ClimaLand.source!(
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

The assumption is that the soil surface temperature
is the same as the temperature at the center of the
first soil layer.
"""
function ClimaLand.component_temperature(
    model::EnergyHydrology{FT},
    Y,
    p,
) where {FT}
    return ClimaLand.Domains.top_center_to_surface(p.soil.T)
end

"""
    ClimaLand.surface_emissivity(
        model::EnergyHydrology{FT},
        Y,
        p,
    ) where {FT}

Returns the surface emissivity field of the
`EnergyHydrology` soil model.
"""
function ClimaLand.surface_emissivity(
    model::EnergyHydrology{FT},
    Y,
    p,
) where {FT}
    return model.parameters.emissivity
end

"""
    ClimaLand.component_specific_humidity(model::EnergyHydrology, Y, p)

a helper function which returns the surface specific humidity for the canopy
model.
"""
function ClimaLand.component_specific_humidity(model::EnergyHydrology, Y, p)
    earth_param_set = get_earth_param_set(model)
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    surface_flux_params =
        LP.surface_fluxes_parameters(model.parameters.earth_param_set)
    T_sfc = component_temperature(model, Y, p)
    Tf_depressed_sfc =
        ClimaLand.Domains.top_center_to_surface(p.soil.Tf_depressed)
    h_sfc = ClimaLand.surface_height(model, Y, p)
    atmos = model.boundary_conditions.top.atmos
    ρ_sfc = @.lazy(
        ClimaLand.compute_ρ_sfc(
            surface_flux_params,
            p.drivers.T,
            p.drivers.P,
            p.drivers.q,
            atmos.h - h_sfc,
            T_sfc,
        ),
    )
    q_sfc = @. lazy(
        soil_specific_humidity(T_sfc, ρ_sfc, Tf_depressed_sfc, earth_param_set),
    )
    return q_sfc
end

function soil_specific_humidity(
    T_sfc::FT,
    ρ_sfc::FT,
    Tf_depressed_sfc::FT,
    earth_param_set,
) where {FT}
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    # Compute q_soil using ice or liquid as appropriate
    if T_sfc > Tf_depressed_sfc # liquid water evaporation
        q_sfc = Thermodynamics.q_vap_saturation(
            thermo_params,
            T_sfc,
            ρ_sfc,
            Thermodynamics.Liquid(),
        )
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
        T_sfc::FT,
        u_star::FT,
        z_0m::FT,
        z_0b::FT,
        K_sfc::FT,
        K_c::FT,
        β_ice::FT,
        Tf_depressed::FT,
    )::FT where {FT}
        g_h = SurfaceFluxes.heat_conductance(
            param_set,
            ζ,
            u_star,
            inputs,
            z_0m,
            z_0b,
            scheme,
        )
        q_air::FT = inputs.q_tot_int - inputs.q_liq_int - inputs.q_ice_int
        qsat_sfc::FT = inputs.q_vap_sfc_guess
        Ẽ_pot::FT = -inputs.ρ_int * g_h * (q_air - qsat_sfc)
        β = FT(1)
        if q_air < qsat_sfc # water loss to atmosphere, adjust β
            if inputs.T_sfc_guess > Tf_depressed # T_sfc is not updated, so T_sfc = T_sfc_guess
                x = 4 * K_sfc * (1 + Ẽ_pot / (4 * K_c))
                β = x / (Ẽ_pot + x)
            else # sublimation, set evaporation to zero
                β = β_ice
            end
        end
        q = β * qsat_sfc + (1 - β) * q_air # q_vap_sfc_guess is already the saturated value
        return q
    end
    # Closure
    FT = eltype(Y)
    earth_param_set = get_earth_param_set(model)
    thermo_params = LP.thermodynamic_parameters(earth_param_set)
    T_sfc = component_temperature(model, Y, p)
    Tf_depressed_sfc =
        ClimaLand.Domains.top_center_to_surface(p.soil.Tf_depressed)
    (; K_sat, ν, θ_r, hydrology_cm, Ω, γ, γT_ref, earth_param_set) =
        model.parameters
    K_sat_sfc = ClimaLand.Domains.top_center_to_surface(K_sat)
    hydrology_cm_sfc = ClimaLand.Domains.top_center_to_surface(hydrology_cm)
    θ_i_sfc = ClimaLand.Domains.top_center_to_surface(Y.soil.θ_i)
    ν_sfc = ClimaLand.Domains.top_center_to_surface(ν)
    θ_r_sfc = ClimaLand.Domains.top_center_to_surface(θ_r)
    θ_l_sfc = p.soil.sfc_scratch
    ClimaLand.Domains.linear_interpolation_to_surface!(
        θ_l_sfc,
        p.soil.θ_l,
        model.domain.fields.z,
        model.domain.fields.Δz_top,
    )
    @. θ_l_sfc = max(θ_l_sfc, θ_r_sfc + sqrt(eps(FT)))
    @. θ_i_sfc = max(θ_i_sfc, FT(0))
    _ρ_liq = LP.ρ_cloud_liq(earth_param_set)

    K_sfc = p.soil.sfc_scratch
    @. K_sfc =
        impedance_factor(θ_i_sfc / (θ_l_sfc + θ_i_sfc - θ_r_sfc), Ω) *
        viscosity_factor(T_sfc, γ, γT_ref) *
        hydraulic_conductivity(
            hydrology_cm_sfc,
            K_sat_sfc,
            effective_saturation(ν_sfc, θ_l_sfc, θ_r_sfc),
        ) *
        _ρ_liq # as a mass flux
    K_c = @. lazy(
        max(
            hydraulic_conductivity(
                hydrology_cm_sfc,
                K_sat_sfc,
                hydrology_cm_sfc.S_c,
            ),
            sqrt(eps(FT)),
        ) * _ρ_liq,
    ) # as a mass flux
    update_q_vap_sfc_field(K_sfc, K_c, β_ice, Tf_depressed) =
        (args...) -> update_q_vap_sfc_at_a_point(
            args...,
            K_sfc,
            K_c,
            β_ice,
            Tf_depressed,
        )
    return @. lazy(
        update_q_vap_sfc_field(
            K_sfc,
            K_c,
            (θ_i_sfc / ν_sfc)^4,
            Tf_depressed_sfc,
        ),
    ) # β_ice = (θ_i_sfc / ν_sfc)^4
end

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

function turbulent_fluxes!(
    dest,
    atmos::AbstractAtmosphericDrivers,
    model::EnergyHydrology,
    Y,
    p,
    t,
)

    T_sfc = component_temperature(model, Y, p)
    q_sfc = component_specific_humidity(model, Y, p)
    roughness_model = surface_roughness_model(model, Y, p)
    update_T_sfc = get_update_surface_temperature_function(model, Y, p)
    update_q_sfc = get_update_surface_humidity_function(model, Y, p)
    h_sfc = surface_height(model, Y, p)
    displ = surface_displacement_height(model, Y, p)
    update_∂T_sfc∂T = get_∂T_sfc∂T_function(model, Y, p)
    update_∂q_sfc∂T = get_∂q_sfc∂T_function(model, Y, p)
    earth_param_set = get_earth_param_set(model)
    momentum_fluxes = Val(return_momentum_fluxes(atmos))
    Tf_depressed_sfc =
        ClimaLand.Domains.top_center_to_surface(p.soil.Tf_depressed)
    gustiness = SurfaceFluxes.ConstantGustinessSpec(atmos.gustiness)
    dest .=
        soil_turbulent_fluxes_at_a_point.(
            momentum_fluxes, # return_extra_fluxes
            ClimaLand.heaviside.(T_sfc, Tf_depressed_sfc), # is_liquid
            p.drivers.P,
            p.drivers.T,
            p.drivers.q, # q_tot
            p.drivers.u,
            atmos.h,
            T_sfc,
            q_sfc,
            roughness_model,
            update_T_sfc,
            update_q_sfc,
            h_sfc,
            displ,
            update_∂T_sfc∂T,
            update_∂q_sfc∂T,
            gustiness,
            earth_param_set,
        )
    return nothing
end

"""
    soil_turbulent_fluxes_at_a_point(return_extra_fluxes, is_liquid, args...;)

This is a wrapper function that allows us to dispatch on the type of `return_extra_fluxes`
as we compute the soil turbulent fluxes pointwise. This is needed because space for the
extra fluxes is only allocated in the cache when running with a `CoupledAtmosphere`.
The function `soil_compute_turbulent_fluxes_at_a_point` does the actual flux computation.

The `return_extra_fluxes` argument indicates whether to return the following:
- momentum fluxes (`ρτxz`, `ρτyz`)
- buoyancy flux (`buoy_flux`)

The field `is_liquid` indicates if the vapor flux is attributed to liquid water evaporating
or due to ice sublimating.
"""
function soil_turbulent_fluxes_at_a_point(
    return_extra_fluxes::Val{false},
    is_liquid,
    args...,
)
    (lhf, shf, vapor_flux, _, _, _, _, _) =
        ClimaLand.compute_turbulent_fluxes_at_a_point(args...)
    return (;
        lhf,
        shf,
        vapor_flux_liq = vapor_flux * is_liquid,
        vapor_flux_ice = vapor_flux * (1 - is_liquid),
    )
end
function soil_turbulent_fluxes_at_a_point(
    return_extra_fluxes::Val{true},
    is_liquid,
    args...,
)
    (lhf, shf, vapor_flux, ∂lhf∂T, ∂shf∂T, ρτxz, ρτyz, buoyancy_flux) =
        ClimaLand.compute_turbulent_fluxes_at_a_point(args...)
    return (;
        lhf,
        shf,
        vapor_flux_liq = vapor_flux * is_liquid,
        vapor_flux_ice = vapor_flux * (1 - is_liquid),
        ρτxz,
        ρτyz,
        buoyancy_flux,
    )
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
        Y.soil.ϑ_l .+ Y.soil.θ_i .* _ρ_ice ./ _ρ_liq, # this line allocates
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

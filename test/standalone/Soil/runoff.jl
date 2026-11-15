import ClimaComms
ClimaComms.@import_required_backends
import ClimaUtilities
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput
using Dates
import ClimaParams as CP
import ClimaLand.Parameters as LP
using ClimaLand
using ClimaLand.Soil.Runoff
using Test
using ClimaCore, NCDatasets
FT = Float32
@testset "Base runoff functionality, FT = $FT" begin
    runoff = Runoff.NoRunoff()
    @test Runoff.runoff_vars(runoff) == (:infiltration,)
    @test Runoff.runoff_var_domain_names(runoff) == (:surface,)
    @test Runoff.runoff_var_types(runoff, FT) == (FT,)

    precip = 5.0
    @test Runoff.subsurface_runoff_source(runoff) == nothing
    struct Foo{FT} <: Runoff.AbstractSoilSource{FT} end
    srcs = (1, 2, 3)
    @test ClimaLand.Soil.append_source(nothing, srcs) == srcs
    @test ClimaLand.Soil.append_source(Foo{FT}(), srcs) == (srcs..., Foo{FT}())
end

@testset "Richards model, no runoff FT =$FT" begin
    domain = ClimaLand.Domains.SphericalShell(;
        radius = FT(6300e3),
        depth = FT(50.0),
        nelements = (101, 15),
        dz_tuple = FT.((5.0, 0.05)),
    )
    surface_space = domain.space.surface
    subsurface_space = domain.space.subsurface
    vg_α = ClimaCore.Fields.ones(subsurface_space) .* FT(0.2)
    hydrology_cm = map(vg_α) do (α)
        FT = typeof(α)
        ClimaLand.Soil.vanGenuchten{FT}(; α = α, n = α + FT(2))
    end
    θ_r = ClimaCore.Fields.zeros(subsurface_space)
    ν = ClimaCore.Fields.zeros(subsurface_space) .+ FT(0.5)
    K_sat = ClimaCore.Fields.zeros(subsurface_space) .+ FT(1e-6)
    S_s = ClimaCore.Fields.zeros(subsurface_space) .+ FT(1e-3)
    soil_params = ClimaLand.Soil.RichardsParameters(;
        hydrology_cm = hydrology_cm,
        ν = ν,
        K_sat = K_sat,
        S_s = S_s,
        θ_r = θ_r,
    )
    lat = ClimaCore.Fields.coordinate_field(surface_space).lat
    precip_field = @. FT(-1e-6) + FT(5e-7) * sin(lat / FT(90.0 * 2π))
    function precip_function(t; fieldvals = precip_field)
        return fieldvals
    end

    precip = TimeVaryingInput(precip_function)
    atmos = ClimaLand.PrescribedPrecipitation{FT, typeof(precip)}(precip)

    noflux = ClimaLand.Soil.WaterFluxBC((p, t) -> 0.0)
    bc = (;
        top = ClimaLand.Soil.RichardsAtmosDrivenFluxBC(atmos),
        bottom = noflux,
    )
    model = ClimaLand.Soil.RichardsModel{FT}(;
        parameters = soil_params,
        domain = domain,
        boundary_conditions = bc,
        sources = (),
    )
    Y, p, t = initialize(model)
    # set initial conditions
    z = ClimaCore.Fields.coordinate_field(domain.space.subsurface).z
    Y.soil.ϑ_l .= FT(0.6) .- FT(0.3 / 50) .* (z .+ FT(50))
    evaluate!(p.drivers.P_liq, atmos.liquid_precip, FT(0))
    set_initial_cache! = make_set_initial_cache(model)
    set_initial_cache!(p, Y, FT(0))
    @test p.soil.top_bc == p.drivers.P_liq
    @test typeof(bc.top.runoff) <: Runoff.NoRunoff
    @test :infiltration ∈ propertynames(p.soil)

    @test Runoff.get_saturated_height(bc.top.runoff, Y, p) isa Nothing
    @test Runoff.get_subsurface_runoff(bc.top.runoff, Y, p) isa Nothing
    @test Runoff.get_surface_runoff(bc.top.runoff, Y, p) isa Nothing
    @test Runoff.get_soil_fsat(bc.top.runoff, Y, p, model.domain.depth) isa
          Nothing

end

@testset "Richards model, TOPMODEL runoff FT =$FT" begin
    domain = ClimaLand.Domains.SphericalShell(;
        radius = FT(6300e3),
        depth = FT(50.0),
        nelements = (101, 15),
        dz_tuple = FT.((5.0, 0.05)),
    )
    surface_space = domain.space.surface
    subsurface_space = domain.space.subsurface
    f_max = ClimaCore.Fields.ones(surface_space) .* FT(0.5)
    f_over = FT(3.28) # 1/m
    R_sb = FT(1.484e-4 / 1000) # m/s
    runoff_model =
        Runoff.TOPMODELRunoff{FT}(; f_over = f_over, f_max = f_max, R_sb = R_sb)
    @test Runoff.runoff_vars(runoff_model) == (
        :infiltration,
        :is_saturated,
        :R_s,
        :R_ss,
        :R_ess,
        :h∇,
        :subsfc_scratch,
    )
    @test Runoff.runoff_var_domain_names(runoff_model) == (
        :surface,
        :subsurface,
        :surface,
        :surface,
        :surface,
        :surface,
        :subsurface,
    )
    @test Runoff.runoff_var_types(runoff_model, FT) ==
          (FT, FT, FT, FT, FT, FT, FT)

    @test runoff_model.f_over == f_over
    @test runoff_model.f_max == f_max
    @test runoff_model.subsurface_source ==
          Runoff.TOPMODELSubsurfaceRunoff{FT}(R_sb, f_over)

    vg_α = ClimaCore.Fields.ones(subsurface_space) .* FT(0.2)
    hydrology_cm = map(vg_α) do (α)
        FT = typeof(α)
        ClimaLand.Soil.vanGenuchten{FT}(; α = α, n = α + FT(2))
    end
    θ_r = ClimaCore.Fields.zeros(subsurface_space)
    ν = ClimaCore.Fields.zeros(subsurface_space) .+ FT(0.5)
    K_sat = ClimaCore.Fields.zeros(subsurface_space) .+ FT(1e-6)
    S_s = ClimaCore.Fields.zeros(subsurface_space) .+ FT(1e-3)
    soil_params = ClimaLand.Soil.RichardsParameters(;
        hydrology_cm = hydrology_cm,
        ν = ν,
        K_sat = K_sat,
        S_s = S_s,
        θ_r = θ_r,
    )
    lat = ClimaCore.Fields.coordinate_field(surface_space).lat
    precip_field = @. FT(-1e-6) + FT(5e-7) * sin(lat / FT(90.0 * 2π))
    function precip_function(t; fieldvals = precip_field)
        return fieldvals
    end

    precip = TimeVaryingInput(precip_function)
    atmos = ClimaLand.PrescribedPrecipitation{FT, typeof(precip)}(precip)

    noflux = ClimaLand.Soil.WaterFluxBC((p, t) -> 0.0)
    bc = (;
        top = ClimaLand.Soil.RichardsAtmosDrivenFluxBC(atmos, runoff_model),
        bottom = noflux,
    )
    model = ClimaLand.Soil.RichardsModel{FT}(;
        parameters = soil_params,
        domain = domain,
        boundary_conditions = bc,
        sources = (),
    )
    Y, p, t = initialize(model)
    @test :R_s ∈ propertynames(p.soil)
    @test :R_ss ∈ propertynames(p.soil)
    @test :R_ess ∈ propertynames(p.soil)
    @test :h∇ ∈ propertynames(p.soil)
    @test :infiltration ∈ propertynames(p.soil)
    @test :is_saturated ∈ propertynames(p.soil)
    @test :subsfc_scratch ∈ propertynames(p.soil)

    # set initial conditions
    z = ClimaCore.Fields.coordinate_field(domain.space.subsurface).z
    Y.soil.ϑ_l .= FT(0.6) .- FT(0.3 / 50) .* (z .+ FT(50))
    evaluate!(p.drivers.P_liq, atmos.liquid_precip, FT(0))
    set_initial_cache! = make_set_initial_cache(model)
    set_initial_cache!(p, Y, FT(0))
    scratch = ClimaCore.zeros(surface_space)
    ClimaCore.Operators.column_integral_definite!(
        scratch,
        ClimaLand.heaviside.(Y.soil.ϑ_l .- model.parameters.ν) .*
        (Y.soil.ϑ_l .- model.parameters.θ_r) ./
        (model.parameters.ν .- model.parameters.θ_r),
    )
    @test scratch == p.soil.h∇
    @test p.soil.R_ss == Runoff.topmodel_ss_flux.(
        runoff_model.subsurface_source.R_sb,
        runoff_model.f_over,
        model.domain.depth .- p.soil.h∇,
    )
    ic = Runoff.soil_infiltration_capacity(model, Y, p)
    @test ic == ClimaCore.Fields.zeros(surface_space) .- FT(1e-6) #Ksat
    @test p.soil.infiltration == @. Runoff.topmodel_surface_infiltration(
        f_max,
        f_over,
        model.domain.depth - p.soil.h∇,
        ic,
        precip_field,
    )
    @test p.soil.R_s == abs.(precip_field .- p.soil.infiltration)

    @test Runoff.get_saturated_height(runoff_model, Y, p) == p.soil.h∇
    @test Runoff.get_subsurface_runoff(runoff_model, Y, p) == p.soil.R_ss
    @test Runoff.get_surface_runoff(runoff_model, Y, p) == p.soil.R_s
    tmp = p.soil.R_s
    tmp .= Runoff.get_soil_fsat(runoff_model, Y, p, model.domain.depth)
    @test tmp == @. runoff_model.f_max *
             exp(-runoff_model.f_over / 2 * (model.domain.depth - p.soil.h∇))

end

@testset "EnergyHydrology model, TOPMODEL runoff FT =$FT" begin
    domain = ClimaLand.Domains.SphericalShell(;
        radius = FT(6300e3),
        depth = FT(50.0),
        nelements = (101, 15),
        dz_tuple = FT.((5.0, 0.05)),
    )
    toml_dict = LP.create_toml_dict(FT)
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        DateTime(2008),
        DateTime(2009),
        domain.space.surface,
        toml_dict,
        FT;
        max_wind_speed = 25.0,
        use_lowres_forcing = true,
    )
    forcing = (; atmos, radiation)
    model = ClimaLand.Soil.EnergyHydrology{FT}(domain, forcing, toml_dict)
    Y, p, t = initialize(model)
    @test :R_s ∈ propertynames(p.soil)
    @test :R_ss ∈ propertynames(p.soil)
    @test :R_ess ∈ propertynames(p.soil)
    @test :h∇ ∈ propertynames(p.soil)
    @test :infiltration ∈ propertynames(p.soil)
    @test :is_saturated ∈ propertynames(p.soil)
    @test :subsfc_scratch ∈ propertynames(p.soil)

    # set initial conditions
    z = ClimaCore.Fields.coordinate_field(domain.space.subsurface).z
    Y.soil.ϑ_l .= FT(0.6) .- FT(0.3 / 50) .* (z .+ FT(50))
    Y.soil.θ_i .= FT(0)
    T = FT(273.17)
    ρc_s = ClimaLand.Soil.volumetric_heat_capacity.(
        Y.soil.ϑ_l,
        Y.soil.θ_i,
        model.parameters.ρc_ds,
        model.parameters.earth_param_set,
    )
    Y.soil.ρe_int .= ClimaLand.Soil.volumetric_internal_energy.(
        Y.soil.θ_i,
        ρc_s,
        T,
        model.parameters.earth_param_set,
    )
    set_initial_cache! = make_set_initial_cache(model)
    set_initial_cache!(p, Y, FT(0))
    surface_space = domain.space.surface
    scratch = ClimaCore.zeros(surface_space)
    runoff_model = model.boundary_conditions.top.runoff
    ClimaCore.Operators.column_integral_definite!(
        scratch,
        ClimaLand.heaviside.(Y.soil.ϑ_l .- model.parameters.ν) .*
        (Y.soil.ϑ_l .- model.parameters.θ_r) ./
        (model.parameters.ν .- model.parameters.θ_r),
    )
    @test scratch == p.soil.h∇
    @test p.soil.R_ss == Runoff.topmodel_ss_flux.(
        runoff_model.subsurface_source.R_sb,
        runoff_model.f_over,
        model.domain.depth .- p.soil.h∇,
    )
    ic = Runoff.soil_infiltration_capacity(model, Y, p)
    @test p.soil.infiltration == @. Runoff.topmodel_surface_infiltration(
        runoff_model.f_max,
        runoff_model.f_over,
        model.domain.depth - p.soil.h∇,
        ic,
        p.drivers.P_liq,
    )
    @test p.soil.R_s == abs.(p.drivers.P_liq .- p.soil.infiltration)
    scratch2 = ClimaCore.zeros(surface_space)
    ClimaCore.Operators.column_integral_definite!(
        scratch2,
        ClimaLand.heaviside.(Y.soil.ϑ_l .- model.parameters.ν) .*
        (Y.soil.ϑ_l .- model.parameters.θ_r) ./
        (model.parameters.ν .- model.parameters.θ_r) .*
        ClimaLand.Soil.volumetric_internal_energy_liq.(
            p.soil.T,
            model.parameters.earth_param_set,
        ),
    )
    @test p.soil.R_ess == scratch2 .* (p.soil.R_ss ./ max.(p.soil.h∇, eps(FT)))
    @test Runoff.get_saturated_height(runoff_model, Y, p) == p.soil.h∇
    @test Runoff.get_subsurface_runoff(runoff_model, Y, p) == p.soil.R_ss
    @test Runoff.get_surface_runoff(runoff_model, Y, p) == p.soil.R_s
    tmp = p.soil.R_s
    tmp .= Runoff.get_soil_fsat(runoff_model, Y, p, model.domain.depth)
    @test tmp == @. runoff_model.f_max *
             exp(-runoff_model.f_over / 2 * (model.domain.depth - p.soil.h∇))
end

@testset "Richards model, Site level runoff FT =$FT" begin
    domain = ClimaLand.Domains.SphericalShell(;
        radius = FT(6300e3),
        depth = FT(50.0),
        nelements = (101, 15),
        dz_tuple = FT.((5.0, 0.05)),
    )
    surface_space = domain.space.surface
    subsurface_space = domain.space.subsurface
    runoff_model = Runoff.SurfaceRunoff()
    @test Runoff.runoff_vars(runoff_model) ==
          (:is_saturated, :R_s, :infiltration, :subsfc_scratch)
    @test Runoff.runoff_var_domain_names(runoff_model) ==
          (:subsurface, :surface, :surface, :subsurface)
    @test Runoff.runoff_var_types(runoff_model, FT) == (FT, FT, FT, FT)


    @test runoff_model.subsurface_source == nothing

    vg_α = ClimaCore.Fields.ones(subsurface_space) .* FT(0.2)
    hydrology_cm = map(vg_α) do (α)
        FT = typeof(α)
        ClimaLand.Soil.vanGenuchten{FT}(; α = α, n = α + FT(2))
    end
    θ_r = ClimaCore.Fields.zeros(subsurface_space)
    ν = ClimaCore.Fields.zeros(subsurface_space) .+ FT(0.5)
    K_sat = ClimaCore.Fields.zeros(subsurface_space) .+ FT(1e-6)
    S_s = ClimaCore.Fields.zeros(subsurface_space) .+ FT(1e-3)
    soil_params = ClimaLand.Soil.RichardsParameters(;
        hydrology_cm = hydrology_cm,
        ν = ν,
        K_sat = K_sat,
        S_s = S_s,
        θ_r = θ_r,
    )
    lat = ClimaCore.Fields.coordinate_field(surface_space).lat
    precip_field = @. FT(-1e-6) + FT(5e-7) * sin(lat / FT(90.0 * 2π))
    function precip_function(t; fieldvals = precip_field)
        return fieldvals
    end

    precip = TimeVaryingInput(precip_function)
    atmos = ClimaLand.PrescribedPrecipitation{FT, typeof(precip)}(precip)

    noflux = ClimaLand.Soil.WaterFluxBC((p, t) -> 0.0)
    bc = (;
        top = ClimaLand.Soil.RichardsAtmosDrivenFluxBC(atmos, runoff_model),
        bottom = noflux,
    )
    model = ClimaLand.Soil.RichardsModel{FT}(;
        parameters = soil_params,
        domain = domain,
        boundary_conditions = bc,
        sources = (),
    )
    Y, p, t = initialize(model)
    @test :R_s ∈ propertynames(p.soil)
    @test :infiltration ∈ propertynames(p.soil)
    @test :is_saturated ∈ propertynames(p.soil)

    # set initial conditions
    z = ClimaCore.Fields.coordinate_field(domain.space.subsurface).z
    Y.soil.ϑ_l .= FT(0.6) .- FT(0.3 / 50) .* (z .+ FT(50))
    evaluate!(p.drivers.P_liq, atmos.liquid_precip, FT(0))
    set_initial_cache! = make_set_initial_cache(model)
    set_initial_cache!(p, Y, FT(0))
    ic = Runoff.soil_infiltration_capacity(model, Y, p)
    @test ic == ClimaCore.Fields.zeros(surface_space) .- FT(1e-6) #Ksat
    @test p.soil.infiltration == Runoff.surface_infiltration.(
        ic,
        precip_field,
        ClimaLand.Domains.top_center_to_surface(p.soil.is_saturated),
    )
    @test p.soil.R_s == abs.(precip_field .- p.soil.infiltration)
    @test Runoff.get_surface_runoff(runoff_model, Y, p) == p.soil.R_s
end

@testset "Gradual frozen infiltration, pointwise FT = $FT" begin
    toml_dict = LP.create_toml_dict(FT)
    fi = Runoff.GradualFrozenInfiltration(toml_dict)
    @test fi isa Runoff.AbstractFrozenInfiltration
    @test fi.ΔT_thaw == FT(1)
    @test fi.ΔT_frozen == FT(3)
    @test fi.ψ_pwp == FT(-153)
    T0 = FT(273.16)
    ff(T) = Runoff.frozen_fraction(FT(T), T0, fi.ΔT_thaw, fi.ΔT_frozen)
    @test ff(T0 + 5) == 0
    @test ff(T0 + 1) == 0
    @test ff(T0 - 1) ≈ FT(0.5)
    # HTESSEL: about 0.15 at 0 °C (IFS CY41R2 Part IV Eq. 8.70)
    @test ff(T0) ≈ (1 - sin(FT(π) / 4)) / 2
    @test ff(T0 - 3) == 1
    @test ff(T0 - 20) == 1
    Ts = range(T0 - 5, T0 + 3; length = 81)
    @test issorted(ff.(Ts); rev = true)
    @test all(0 .<= ff.(Ts) .<= 1)
    @test eltype(ff.(Ts)) == FT
    # the factor multiplying K_sat blends 1 and K_r(ψ_pwp)
    @test Runoff.gradual_frozen_factor(FT(0), FT(1e-6)) == 1
    @test Runoff.gradual_frozen_factor(FT(1), FT(1e-6)) == FT(1e-6)
    @test Runoff.gradual_frozen_factor(FT(0.5), FT(0)) == FT(0.5)
    # Overriding the global parameters
    fi2 = Runoff.GradualFrozenInfiltration(toml_dict; ΔT_thaw = FT(0.5))
    @test fi2.ΔT_thaw == FT(0.5)
    @test fi2.ΔT_frozen == FT(3)
    # Default option of TOPMODELRunoff is the ice impedance
    f_max = ClimaCore.Fields.zeros(
        ClimaLand.Domains.Column(; zlim = FT.((-1.0, 0.0)), nelements = 5).space.surface,
    )
    r0 = Runoff.TOPMODELRunoff(toml_dict; f_max)
    @test r0.frozen_infiltration isa Runoff.IceImpedanceInfiltration
    @test Runoff.frozen_infiltration_model(r0) isa
          Runoff.IceImpedanceInfiltration
    r1 = Runoff.TOPMODELRunoff(toml_dict; f_max, frozen_infiltration = fi)
    @test Runoff.frozen_infiltration_model(r1) === fi
    r2 = Runoff.TOPMODELRunoff{FT}(;
        f_over = FT(3.28),
        f_max,
        R_sb = FT(1e-4),
        frozen_infiltration = fi,
    )
    @test r2.frozen_infiltration === fi
    @test Runoff.frozen_infiltration_model(Runoff.SurfaceRunoff()) isa
          Runoff.IceImpedanceInfiltration
end

@testset "EnergyHydrology model, TOPMODEL runoff with gradual frozen infiltration FT =$FT" begin
    domain = ClimaLand.Domains.SphericalShell(;
        radius = FT(6300e3),
        depth = FT(50.0),
        nelements = (10, 15),
        dz_tuple = FT.((5.0, 0.05)),
    )
    toml_dict = LP.create_toml_dict(FT)
    atmos, radiation = ClimaLand.prescribed_forcing_era5(
        DateTime(2008),
        DateTime(2009),
        domain.space.surface,
        toml_dict,
        FT;
        max_wind_speed = 25.0,
        use_lowres_forcing = true,
    )
    forcing = (; atmos, radiation)
    fi = Runoff.GradualFrozenInfiltration(toml_dict)
    f_max = ClimaLand.Soil.topmodel_fmax(domain.space.surface, FT)
    runoff = Runoff.TOPMODELRunoff(toml_dict; f_max, frozen_infiltration = fi)
    model_default =
        ClimaLand.Soil.EnergyHydrology{FT}(domain, forcing, toml_dict)
    model =
        ClimaLand.Soil.EnergyHydrology{FT}(domain, forcing, toml_dict; runoff)
    @test model.boundary_conditions.top.runoff.frozen_infiltration === fi

    T0 = LP.T_freeze(model.parameters.earth_param_set)
    surface_space = domain.space.surface
    lat = ClimaCore.Fields.coordinate_field(domain.space.subsurface).lat
    function init(model)
        Y, p, t = initialize(model)
        z = ClimaCore.Fields.coordinate_field(domain.space.subsurface).z
        Y.soil.ϑ_l .= FT(0.6) .- FT(0.3 / 50) .* (z .+ FT(50))
        # Ice everywhere, and temperatures from T0 - 6 K to T0 + 6 K, so
        # that the default ice impedance is large everywhere while the
        # frozen fraction spans [0, 1].
        Y.soil.θ_i .= FT(0.1)
        T = @. T0 + FT(6) * sind(lat)
        ρc_s = ClimaLand.Soil.volumetric_heat_capacity.(
            Y.soil.ϑ_l,
            Y.soil.θ_i,
            model.parameters.ρc_ds,
            model.parameters.earth_param_set,
        )
        Y.soil.ρe_int .= ClimaLand.Soil.volumetric_internal_energy.(
            Y.soil.θ_i,
            ρc_s,
            T,
            model.parameters.earth_param_set,
        )
        set_initial_cache! = make_set_initial_cache(model)
        set_initial_cache!(p, Y, FT(0))
        return Y, p
    end
    # Note: models built on the same domain share cache storage, so the
    # default model is initialized only after the checks on `model`.
    Y, p = init(model)

    # The 3-argument method uses the option of the runoff model
    ic = similar(p.soil.R_s)
    ic .= Runoff.soil_infiltration_capacity(model, Y, p)
    @test ic == Runoff.soil_infiltration_capacity(fi, model, Y, p)
    # Ice impedance capacity of the same state
    ic_default = similar(ic)
    ic_default .= Runoff.soil_infiltration_capacity(
        Runoff.IceImpedanceInfiltration(),
        model,
        Y,
        p,
    )
    # The infiltration uses the gradual infiltration capacity
    Runoff.update_infiltration_water_flux!(
        p,
        runoff,
        p.drivers.P_liq,
        Y,
        FT(0),
        model,
    )
    # (the surface runoff uses the water table height including ice; p.soil.h∇
    # holds the liquid-only one after the update)
    h∇_ice = ClimaCore.zeros(surface_space)
    (; ν, θ_r) = model.parameters
    ClimaCore.Operators.column_integral_definite!(
        h∇_ice,
        ClimaLand.heaviside.(Y.soil.ϑ_l .+ Y.soil.θ_i .- ν) .*
        (Y.soil.ϑ_l .+ Y.soil.θ_i .- θ_r) ./ (ν .- θ_r),
    )
    infil = similar(ic)
    infil .= Runoff.topmodel_surface_infiltration.(
        runoff.f_max,
        runoff.f_over,
        model.domain.depth .- h∇_ice,
        ic,
        p.drivers.P_liq,
    )
    @test Array(parent(p.soil.infiltration)) ≈ Array(parent(infil))
    # Expected values on the top layer
    (; K_sat, hydrology_cm, γ, γT_ref) = model.parameters
    top(x) = ClimaLand.Domains.top_center_to_surface(x)
    T_sfc = similar(ic)
    T_sfc .= top(p.soil.T)
    ff = similar(ic)
    ff .= Runoff.frozen_fraction.(T_sfc, T0, fi.ΔT_thaw, fi.ΔT_frozen)
    K_r_pwp = similar(p.soil.T)
    K_r_pwp .= ClimaLand.Soil.hydraulic_conductivity.(
        hydrology_cm,
        FT(1),
        ClimaLand.Soil.inverse_matric_potential.(hydrology_cm, fi.ψ_pwp),
    )
    K_max = similar(p.soil.T)
    K_max .= .-K_sat .* ClimaLand.Soil.viscosity_factor.(p.soil.T, γ, γT_ref)
    expected = similar(ic)
    expected .= top(K_max) .* Runoff.gradual_frozen_factor.(ff, top(K_r_pwp))
    @test Array(parent(ic)) ≈ Array(parent(expected))
    # Unfrozen top layer (T ≥ T0 + 1 K): no reduction despite the ice, and
    # larger than with the ice impedance. Fully frozen (T ≤ T0 - 3 K):
    # capacity reduced to the conductivity at the wilting point.
    ff_a = Array(parent(ff))
    ic_a = Array(parent(ic))
    icd_a = Array(parent(ic_default))
    Kmax_a = Array(parent(top(K_max)))
    Krp_a = Array(parent(top(K_r_pwp)))
    warm = ff_a .== 0
    cold = ff_a .== 1
    @test any(warm) && any(cold) && any(0 .< ff_a .< 1)
    @test ic_a[warm] ≈ Kmax_a[warm]
    @test all(ic_a[warm] .< icd_a[warm])
    @test ic_a[cold] ≈ Kmax_a[cold] .* Krp_a[cold]
    @test all(0 .<= Krp_a .< FT(1e-3))
    # The interior conductivity does not depend on the option
    K_gradual = copy(Array(parent(p.soil.K)))
    T_gradual = copy(Array(parent(p.soil.T)))
    Yd, pd = init(model_default)
    @test Array(parent(pd.soil.K)) == K_gradual
    @test Array(parent(pd.soil.T)) == T_gradual
    # The default model keeps the ice impedance
    @test model_default.boundary_conditions.top.runoff.frozen_infiltration isa
          Runoff.IceImpedanceInfiltration
    ic_d = similar(ic)
    ic_d .= Runoff.soil_infiltration_capacity(model_default, Yd, pd)
    @test ic_d == ic_default
end

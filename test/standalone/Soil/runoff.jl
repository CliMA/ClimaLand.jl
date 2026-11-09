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

    f_frz = Runoff.frozen_impermeable_fraction
    @test f_frz(FT(0)) ≈ 0 atol = eps(FT)
    @test f_frz(FT(1)) ≈ 1
    @test f_frz(FT(-0.1)) == f_frz(FT(0))
    @test f_frz(FT(1.2)) == f_frz(FT(1))
    @test issorted(f_frz.(FT.(0:0.1:1)))
    @test f_frz(FT(0.5)) < FT(0.2)

    @test Runoff.is_frozen(FT(0.3), FT(0.1), FT(0.05))
    @test !Runoff.is_frozen(FT(0.1), FT(0.3), FT(0.05))
    @test !Runoff.is_frozen(FT(0), FT(0.3), FT(0.05))
    R_sb, f_over = FT(1.484e-7), FT(3.28)
    @test Runoff.perched_drainage(R_sb, f_over, FT(1), FT(0)) == 0
    @test Runoff.perched_drainage(R_sb, f_over, FT(1), FT(0.5)) ==
          Runoff.topmodel_ss_flux(R_sb, f_over, FT(0.5))
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
        :ss_sink,
        :z∇,
    )
    @test Runoff.runoff_var_domain_names(runoff_model) == (
        :surface,
        :subsurface,
        :surface,
        :surface,
        :surface,
        :surface,
        :subsurface,
        :subsurface,
        :surface,
    )
    @test Runoff.runoff_var_types(runoff_model, FT) ==
          (FT, FT, FT, FT, FT, FT, FT, FT, FT)

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
    @test Array(parent(p.soil.R_ss)) ≈ Array(
        parent(
            Runoff.topmodel_ss_flux.(
                runoff_model.subsurface_source.R_sb,
                runoff_model.f_over,
                model.domain.depth .- p.soil.h∇,
            ),
        ),
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
    @test Array(parent(p.soil.z∇)) ≈
          Array(parent(max.(model.domain.depth .- p.soil.h∇, 0)))
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
    @test Array(parent(p.soil.R_ss)) ≈ Array(
        parent(
            Runoff.topmodel_ss_flux.(
                runoff_model.subsurface_source.R_sb,
                runoff_model.f_over,
                model.domain.depth .- p.soil.h∇,
            ),
        ),
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
    @test Array(parent(p.soil.R_ess)) ≈
          Array(parent(scratch2 .* (p.soil.R_ss ./ max.(p.soil.h∇, eps(FT)))))
    @test Runoff.get_saturated_height(runoff_model, Y, p) == p.soil.h∇
    @test Array(parent(p.soil.z∇)) ≈
          Array(parent(max.(model.domain.depth .- p.soil.h∇, 0)))
    @test Runoff.get_subsurface_runoff(runoff_model, Y, p) == p.soil.R_ss
    @test Runoff.get_surface_runoff(runoff_model, Y, p) == p.soil.R_s
    tmp = p.soil.R_s
    tmp .= Runoff.get_soil_fsat(runoff_model, Y, p, model.domain.depth)
    @test tmp == @. runoff_model.f_max *
             exp(-runoff_model.f_over / 2 * (model.domain.depth - p.soil.h∇))

    # Water perched on frozen layers drains from the frost table
    ν = model.parameters.ν
    θ_r = model.parameters.θ_r
    frozen = z .< FT(-2)
    @. Y.soil.θ_i = ifelse(frozen, FT(0.8) * (ν - θ_r), FT(0))
    @. Y.soil.ϑ_l = ifelse(frozen, θ_r + FT(0.01), ν + FT(0.01))
    @. Y.soil.ρe_int = ClimaLand.Soil.volumetric_internal_energy(
        Y.soil.θ_i,
        ClimaLand.Soil.volumetric_heat_capacity(
            Y.soil.ϑ_l,
            Y.soil.θ_i,
            model.parameters.ρc_ds,
            model.parameters.earth_param_set,
        ),
        ifelse(frozen, FT(265), FT(275)),
        model.parameters.earth_param_set,
    )
    set_initial_cache!(p, Y, FT(0))
    frozen_layers = Array(parent(frozen)) .== 1
    @test all(
        Runoff.is_frozen.(
            Array(parent(Y.soil.θ_i))[frozen_layers],
            Array(parent(Y.soil.ϑ_l))[frozen_layers],
            Array(parent(θ_r))[frozen_layers],
        ),
    )
    integral = ClimaCore.zeros(surface_space)
    ClimaCore.Operators.column_integral_definite!(integral, p.soil.ss_sink)
    @test Array(parent(integral)) ≈ Array(parent(p.soil.R_ss))
    @test all(iszero, Array(parent(p.soil.ss_sink))[frozen_layers])
    deep_table = Runoff.topmodel_ss_flux.(
        runoff_model.subsurface_source.R_sb,
        runoff_model.f_over,
        model.domain.depth .- p.soil.h∇,
    )
    @test all(Array(parent(p.soil.R_ss)) .> 100 .* Array(parent(deep_table)))
    # The water table is that of the perched water: at the surface here
    @test all(iszero, Array(parent(p.soil.z∇)))
    @test Runoff.get_water_table_depth(runoff_model, Y, p) === p.soil.z∇

    # A saturated surface layer takes no water: all liquid input runs off
    ClimaLand.Domains.top_center_to_surface(Y.soil.ϑ_l) .=
        ClimaLand.Domains.top_center_to_surface(model.parameters.ν) .+ FT(0.05)
    set_initial_cache!(p, Y, FT(0))
    @test all(iszero, parent(p.soil.infiltration))
    @test p.soil.R_s == abs.(p.drivers.P_liq)
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

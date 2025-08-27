using Test
import ClimaComms
ClimaComms.@import_required_backends
using Statistics
using ClimaCore
import ClimaParams
using ClimaLand
using ClimaLand.Soil
using ClimaLand.Domains: HybridBox, SphericalShell, Column
import ClimaLand.Parameters as LP
import ClimaParams

for FT in (Float32, Float64)
    default_params_filepath =
        joinpath(pkgdir(ClimaLand), "toml", "default_parameters.toml")
    toml_dict = LP.create_toml_dict(FT, default_params_filepath)
    @testset "WVector usage in gradient, FT = $FT" begin
        # In our tendency, we want to set a boundary condition on flux, which is a gradient
        # of a scalar. This should be a covariant vector. In the case of no topography,
        # the Wvector points in the same direction as the covariant vector.

        # However, they are scaled differently. Check here that supplying a WVector
        # to Operators.DivergenceF2C(top = Operators.SetValue...,) works - i.e.
        # that it is converted internally to a Covariant Vector, and we get the
        # correct output.

        domain = SphericalShell(;
            radius = FT(1.0),
            depth = FT(1.0),
            nelements = (1, 2),
        )

        coords = ClimaLand.Domains.coordinates(domain) # center coordinates
        ϕ = (FT(90) .- coords.subsurface.lat)
        θ = coords.subsurface.long
        r = coords.subsurface.z .+ FT(1.0)

        #=
        # For my own sanity, convert lat/lon to usual spherical coords, to cartesian components,
        # and make sure that agrees with CC.
        scalar_x = @.  r*cos(θ*π/180)*sin(ϕ*π/180)
        scalar_y = @.  r*sin(θ*π/180)*sin(ϕ*π/180)
        scalar_z = @. r* cos(ϕ*π/180)
        geom = ClimaCore.Geometry.SphericalGlobalGeometry{FT}(FT(1.0))
        #\vec{r} = r r^
        cartesian = @. ClimaCore.Geometry.CartesianVector(ClimaCore.Geometry.UVWVector(0.0,0.0, r))
        cartesian_points = ClimaCore.Geometry.CartesianPoint.(coords, Ref(geom))
        =#

        scalar = @. FT(r^3 / 1e6)
        Δr = FT(0.5)
        upper_bc = ClimaCore.Geometry.Covariant3Vector(FT(2.0 * Δr)) # r r^
        lower_bc = ClimaCore.Geometry.Covariant3Vector(FT(1.0 * Δr)) # r r^
        gradc2f = ClimaCore.Operators.GradientC2F(
            top = ClimaCore.Operators.SetGradient(upper_bc),
            bottom = ClimaCore.Operators.SetGradient(lower_bc),
        )
        foo = @. gradc2f(scalar)
        divf2c = ClimaCore.Operators.DivergenceF2C()

        div_cov = @. (divf2c(foo))


        upper_bc = ClimaCore.Geometry.Contravariant3Vector(FT(2.0 * Δr)) # r r^
        lower_bc = ClimaCore.Geometry.Contravariant3Vector(FT(1.0 * Δr)) # r r^
        gradc2f = ClimaCore.Operators.GradientC2F(
            top = ClimaCore.Operators.SetGradient(upper_bc),
            bottom = ClimaCore.Operators.SetGradient(lower_bc),
        )
        foo = @. gradc2f(scalar)
        divf2c = ClimaCore.Operators.DivergenceF2C()

        div_con = @. (divf2c(foo))

        upper_bc = ClimaCore.Geometry.WVector(FT(2.0)) # r r^
        lower_bc = ClimaCore.Geometry.WVector(FT(1.0)) # r r^
        gradc2f = ClimaCore.Operators.GradientC2F(
            top = ClimaCore.Operators.SetGradient(upper_bc),
            bottom = ClimaCore.Operators.SetGradient(lower_bc),
        )
        foo = @. gradc2f(scalar)
        divf2c = ClimaCore.Operators.DivergenceF2C()

        div_w = @. (divf2c(foo))
        @test parent(div_w) == parent(div_cov)
        @test all(abs.(parent(div_w) .- parent(div_cov)) .< eps(FT))

        # Unclear why this does not work!
        # upper_bc =ClimaCore.Geometry.Covariant3Vector(2.0*Δr) # r r^
        # lower_bc =ClimaCore.Geometry.Covariant3Vector(1.0*Δr) # r r^
        # gradc2f = ClimaCore.Operators.GradientC2F()
        # divf2c = ClimaCore.Operators.DivergenceF2C(
        #     top = ClimaCore.Operators.SetValue(upper_bc),
        #     bottom = ClimaCore.Operators.SetValue(lower_bc),
        # )

        #div_rpt = unique(parent(@. (divf2c(gradc2f(scalar)))))
    end

    @testset "Test RRE state to flux BC calculations, FT = $FT" begin
        ν = FT(0.495)
        K_sat = FT(0.0443 / 3600 / 100) # m/s
        S_s = FT(1e-3) #inverse meters
        vg_n = FT(2.0)
        vg_α = FT(2.6) # inverse meters
        hcm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
        θ_r = FT(0)
        zmax = FT(0)
        zmin = FT(-10)
        nelems = 50
        soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
        z = ClimaCore.Fields.coordinate_field(soil_domain.space.subsurface).z
        Δz = FT(abs(zmax - zmin) / nelems / 2.0)

        top_Δz, bottom_Δz = ClimaLand.Domains.get_Δz(z)
        @test (mean(abs.(parent(top_Δz) .- Δz)) < eps(FT))
        @test (mean(abs.(parent(bottom_Δz) .- Δz)) < 2 * eps(FT))

        ϑ_bc = FT(ν / 2)
        ϑ_c = FT(ν / 3)

        K_c = Soil.hydraulic_conductivity(
            hcm,
            K_sat,
            Soil.effective_saturation(ν, ϑ_c, θ_r),
        )

        ψ_bc = pressure_head(hcm, θ_r, ϑ_bc, ν, S_s)
        ψ_c = pressure_head(hcm, θ_r, ϑ_c, ν, S_s)

        flux_int = diffusive_flux(K_c, ψ_bc + Δz, ψ_c, Δz)
        flux_expected = -K_c * ((ψ_bc - ψ_c + Δz) / Δz)

        @test abs(flux_expected - flux_int) < eps(FT)
    end

    @testset "Test heat state to flux BC calculations, FT = $FT" begin
        ν = FT(0.495)
        K_sat = FT(0.0443 / 3600 / 100) # m/s
        S_s = FT(1e-3) #inverse meters
        vg_n = FT(2.0)
        vg_α = FT(2.6) # inverse meters
        hydrology_cm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
        θ_r = FT(0.1)
        ν_ss_om = FT(0.0)
        ν_ss_quartz = FT(1.0)
        ν_ss_gravel = FT(0.0)

        parameters = Soil.EnergyHydrologyParameters(
            FT;
            ν,
            ν_ss_om,
            ν_ss_quartz,
            ν_ss_gravel,
            hydrology_cm,
            K_sat,
            S_s,
            θ_r,
        )

        zmax = FT(0)
        zmin = FT(-1)
        nelems = 200
        Δz = abs(zmax - zmin) / nelems

        ϑ_bc = ν / 2
        ϑ_c = ν / 3
        θ_i = FT(0)
        T_bc = FT(298)
        T_c = FT(290)

        κ_c =
            thermal_conductivity.(
                parameters.κ_dry,
                kersten_number.(
                    θ_i,
                    relative_saturation.(ϑ_c, θ_i, ν),
                    parameters.α,
                    parameters.β,
                    ν_ss_om,
                    ν_ss_quartz,
                    ν_ss_gravel,
                ),
                κ_sat.(
                    ϑ_c,
                    θ_i,
                    parameters.κ_sat_unfrozen,
                    parameters.κ_sat_frozen,
                ),
            )

        flux_int = diffusive_flux(κ_c, T_bc, T_c, Δz)
        flux_expected = -κ_c * (T_bc - T_c) / Δz

        @test abs(flux_expected - flux_int) < eps(FT)
    end

    @testset "Test water+energy free drainage BC, FT = $FT" begin
        ν = FT(0.495)
        K_sat = FT(0.0443 / 3600 / 100) # m/s
        S_s = FT(1e-3) #inverse meters
        vg_n = FT(2.0)
        vg_α = FT(2.6) # inverse meters
        hydrology_cm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
        θ_r = FT(0.1)
        ν_ss_om = FT(0.0)
        ν_ss_quartz = FT(1.0)
        ν_ss_gravel = FT(0.0)

        parameters = Soil.EnergyHydrologyParameters(
            FT;
            ν,
            ν_ss_om,
            ν_ss_quartz,
            ν_ss_gravel,
            hydrology_cm,
            K_sat,
            S_s,
            θ_r,
        )

        zmax = FT(0)
        zmin = FT(-1)
        nelems = 20
        domain = Column(; zlim = (zmin, zmax), nelements = nelems)

        zero_water_flux_bc = WaterFluxBC((p, t) -> 0.0)
        zero_heat_flux_bc = HeatFluxBC((p, t) -> 0.0)
        boundary_conditions = (;
            top = WaterHeatBC(;
                water = zero_water_flux_bc,
                heat = zero_heat_flux_bc,
            ),
            bottom = EnergyWaterFreeDrainage(),
        )
        sources = ()
        soil = Soil.EnergyHydrology{FT}(;
            parameters,
            domain,
            boundary_conditions,
            sources,
        )

        Y, p, coords = initialize(soil)
        Y.soil.ϑ_l .= ν / 2
        Y.soil.θ_i .= 0
        T = coords.subsurface.z .* 10 .+ FT(290.0)
        ρc_s =
            Soil.volumetric_heat_capacity.(
                Y.soil.ϑ_l,
                Y.soil.θ_i,
                parameters.ρc_ds,
                parameters.earth_param_set,
            )
        Y.soil.ρe_int .=
            Soil.volumetric_internal_energy.(
                Y.soil.θ_i,
                ρc_s,
                T,
                parameters.earth_param_set,
            )
        set_initial_cache! = ClimaLand.make_set_initial_cache(soil)
        set_initial_cache!(p, Y, 0.0)
        flux_expected_water = -1 .* ClimaCore.Fields.level(p.soil.K, 1)
        T1 = ClimaCore.Fields.level(p.soil.T, 1)
        flux_expected_energy =
            flux_expected_water .*
            Soil.volumetric_internal_energy_liq.(T1, parameters.earth_param_set)
        @test all(parent(p.soil.bottom_bc.heat) .≈ parent(flux_expected_energy))
        @test all(parent(p.soil.bottom_bc.water) .≈ parent(flux_expected_water))

    end
end


@testset "Boundary vars" begin
    # AtmosDrivenFluxBC method tested in `test/standalone/Soil/climate_drivers.jl`
    # Currently no other methods exist besides the default, which we test here
    FT = Float32
    default_params_filepath =
        joinpath(pkgdir(ClimaLand), "toml", "default_parameters.toml")
    toml_dict = LP.create_toml_dict(FT, default_params_filepath)
    ν = FT(0.495)
    K_sat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(2.0)
    vg_α = FT(2.6) # inverse meters
    hydrology_cm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
    θ_r = FT(0.1)
    ν_ss_om = FT(0.0)
    ν_ss_quartz = FT(1.0)
    ν_ss_gravel = FT(0.0)
    κ_minerals = FT(2.5)
    ρp = FT(2.66 / 1e3 * 1e6)
    ρc_ds = @. FT(2e6 * (1.0 - ν))
    parameters = Soil.EnergyHydrologyParameters(
        FT;
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm,
        K_sat,
        S_s,
        θ_r,
    )
    zmax = FT(0)
    zmin = FT(-1)
    nelems = 200
    Δz = abs(zmax - zmin) / nelems
    soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
    top_heat_flux_bc = HeatFluxBC((p, t) -> -1.0e-8)
    top_water_flux_bc = WaterFluxBC((p, t) -> -1.0e-8)
    bot_heat_flux_bc = HeatFluxBC((p, t) -> -2.0e-8)
    bot_water_flux_bc = WaterFluxBC((p, t) -> -2.0e-8)
    sources = ()
    # Use the same BCs for RRE and heat
    boundary_fluxes = (;
        top = WaterHeatBC(; water = top_water_flux_bc, heat = top_heat_flux_bc),
        bottom = WaterHeatBC(;
            water = bot_water_flux_bc,
            heat = bot_heat_flux_bc,
        ),
    )

    parameters = Soil.EnergyHydrologyParameters(
        FT;
        ν,
        ν_ss_om,
        ν_ss_quartz,
        ν_ss_gravel,
        hydrology_cm,
        K_sat = FT(0),
        S_s,
        θ_r,
    )
    energy_hydrology = Soil.EnergyHydrology{FT}(;
        parameters = parameters,
        domain = soil_domain,
        boundary_conditions = boundary_fluxes,
        sources = sources,
    )

    @test Soil.boundary_vars(boundary_fluxes.top, ClimaLand.TopBoundary()) ==
          (:top_bc, :top_bc_wvec)
    @test Soil.boundary_var_domain_names(
        boundary_fluxes.top,
        ClimaLand.TopBoundary(),
    ) == (:surface, :surface)
    @test Soil.boundary_var_types(
        energy_hydrology,
        boundary_fluxes.top,
        ClimaLand.TopBoundary(),
    ) == (
        NamedTuple{(:water, :heat), Tuple{Float32, Float32}},
        ClimaCore.Geometry.WVector{Float32},
    )

    @test Soil.boundary_vars(
        boundary_fluxes.bottom,
        ClimaLand.BottomBoundary(),
    ) == (:bottom_bc, :bottom_bc_wvec)
    @test Soil.boundary_var_domain_names(
        boundary_fluxes.bottom,
        ClimaLand.BottomBoundary(),
    ) == (:surface, :surface)
    @test Soil.boundary_var_types(
        energy_hydrology,
        boundary_fluxes.bottom,
        ClimaLand.BottomBoundary(),
    ) == (
        NamedTuple{(:water, :heat), Tuple{Float32, Float32}},
        ClimaCore.Geometry.WVector{Float32},
    )
    Y, p, cds = initialize(energy_hydrology)
    Y.soil.ϑ_l .= ν
    Y.soil.θ_i .= 0
    Y.soil.ρe_int .= 3e7
    update_boundary_vars! = make_update_boundary_fluxes(energy_hydrology)
    update_boundary_vars!(p, Y, FT(0))
    f = similar(p.soil.top_bc.water)
    fill!(ClimaCore.Fields.field_values(f), FT(top_heat_flux_bc.bc(p, FT(0))))
    @test p.soil.top_bc.water == f
    @test p.soil.top_bc.heat == f

    f = similar(p.soil.bottom_bc.water)
    fill!(ClimaCore.Fields.field_values(f), FT(bot_heat_flux_bc.bc(p, FT(0))))
    @test p.soil.bottom_bc.water == f
    @test p.soil.bottom_bc.heat == f


    boundary_fluxes = (; top = top_water_flux_bc, bottom = bot_water_flux_bc)

    parameters = Soil.RichardsParameters(;
        ν = ν,
        hydrology_cm,
        K_sat = FT(0),
        S_s = S_s,
        θ_r = θ_r,
    )
    rre = Soil.RichardsModel{FT}(;
        parameters = parameters,
        domain = soil_domain,
        boundary_conditions = boundary_fluxes,
        sources = sources,
    )

    @test Soil.boundary_vars(boundary_fluxes.top, ClimaLand.TopBoundary()) ==
          (:top_bc, :top_bc_wvec)
    @test Soil.boundary_var_domain_names(
        boundary_fluxes.top,
        ClimaLand.TopBoundary(),
    ) == (:surface, :surface)
    @test Soil.boundary_var_types(
        rre,
        boundary_fluxes.top,
        ClimaLand.TopBoundary(),
    ) == (FT, ClimaCore.Geometry.WVector{Float32})

    @test Soil.boundary_vars(boundary_fluxes.top, ClimaLand.BottomBoundary()) ==
          (:bottom_bc, :bottom_bc_wvec)
    @test Soil.boundary_var_domain_names(
        boundary_fluxes.top,
        ClimaLand.BottomBoundary(),
    ) == (:surface, :surface)
    @test Soil.boundary_var_types(
        rre,
        boundary_fluxes.bottom,
        ClimaLand.BottomBoundary(),
    ) == (FT, ClimaCore.Geometry.WVector{FT})

    Y, p, cds = initialize(rre)
    Y.soil.ϑ_l .= ν
    update_boundary_vars! = make_update_boundary_fluxes(rre)
    update_boundary_vars!(p, Y, FT(0))
    f = similar(p.soil.top_bc)
    fill!(
        ClimaCore.Fields.field_values(f),
        FT(boundary_fluxes.top.bc(p, FT(0))),
    )
    @test p.soil.top_bc == f
    f = similar(p.soil.bottom_bc)
    fill!(
        ClimaCore.Fields.field_values(f),
        FT(boundary_fluxes.bottom.bc(p, FT(0))),
    )
    @test p.soil.bottom_bc == f
end

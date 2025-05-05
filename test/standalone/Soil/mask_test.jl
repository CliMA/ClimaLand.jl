import ClimaCore
using ClimaLand
using ClimaLand.Soil
using Test

if pkgversion(ClimaCore) >= v"0.14.30"
    FT = Float64
    domain = ClimaLand.ModelSetup.global_domain(FT)

    # Soil model setup
    ν = FT(0.495)
    K_sat = FT(0.0443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(2.0)
    vg_α = FT(2.6) # inverse meters
    hcm = vanGenuchten{FT}(; α = vg_α, n = vg_n)
    θ_r = FT(0)
    top_flux_bc = WaterFluxBC((p, t) -> 0.0)
    bot_flux_bc = WaterFluxBC((p, t) -> 0.0)
    sources = ()
    boundary_fluxes = (; top = top_flux_bc, bottom = bot_flux_bc)
    params = Soil.RichardsParameters(;
        ν = ν,
        hydrology_cm = hcm,
        K_sat = K_sat,
        S_s = S_s,
        θ_r = θ_r,
    )

    soil = Soil.RichardsModel{FT}(;
        parameters = params,
        domain = domain,
        boundary_conditions = boundary_fluxes,
        sources = sources,
    )

    Y, p, coords = initialize(soil) # Y is the state vector, initially all zero
    # test 2d
    Y .= 2.0 # set to something nonzero using FieldVector broadcasting, which operates on masked areas too.
    sfc_field = ClimaLand.top_center_to_surface(Y.soil.ϑ_l)
    sfc_field .= 1.0 # Field broadcasting respects masks
    @test extrema(sfc_field) == (1.0, 2.0)

    # test 3d
    Y.soil.ϑ_l .= 1.0
    @test extrema(Y.soil.ϑ_l) == (1.0, 2.0)


    # More complex function
    t0 = FT(0)
    set_initial_cache! = make_set_initial_cache(soil)
    set_initial_cache!(p, Y, t0)
    dY = similar(Y)
    @. dY = 2.0

    imp_tendency! = make_imp_tendency(soil)
    imp_tendency!(dY, Y, p, t0)

    binary_mask = .~parent(domain.space.surface.grid.mask.is_active)[:]
    # Test that the masked parts of dY did not update and are still equal to 2.0
    @test extrema(parent(dY.soil.ϑ_l)[:, 1, 1, 1, binary_mask]) == (2.0, 2.0)
end

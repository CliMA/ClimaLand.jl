using ClimaLSM
using Test
using ClimaCore
import CLIMAParameters as CP
using ClimaLSM.Domains: HybridBox, SphericalShell
using ClimaLSM.Soil
import ClimaLSM
import ClimaLSM.Parameters as LSMP
using ClimaLSM.Regridder: regrid_netcdf_to_field
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))


@testset "Base runoff functionality" begin
    runoff = ClimaLSM.Soil.NoRunoff()
    precip = 5.0
    @test ClimaLSM.Soil.soil_surface_infiltration(runoff, precip) == precip
    @test ClimaLSM.Soil.subsurface_runoff_source(runoff) == nothing
    struct Foo{Float64} <: ClimaLSM.Soil.AbstractSoilSource{Float64} end
    srcs = (1, 2, 3)
    @test ClimaLSM.Soil.append_source(nothing, srcs) == srcs
    @test ClimaLSM.Soil.append_source(Foo{Float64}(), srcs) ==
          (srcs..., Foo{Float64}())
end


@test "TOPMODEL runoff" begin
    FT = Float64
    earth_param_set = create_lsm_parameters(FT)
    r_earth = LSMP.r_earth(earth_param_set)
    ν = FT(0.495)
    K_sat = FT(0.00443 / 3600 / 100) # m/s
    S_s = FT(1e-3) #inverse meters
    vg_n = FT(2.0)
    vg_α = FT(2.6) # inverse meters
    vg_m = 1 - 1 / vg_n
    hcm = vanGenuchten(; α = vg_α, n = vg_n)
    θ_r = FT(0)
    soil_domain = SphericalShell(;
        radius = r_earth,
        depth = FT(5.0),
        nelements = (1, 10),
        npolynomial = 1,
    )
    bot_flux_bc = FluxBC((p, t) -> eltype(t)(0.0))
    precip = (t) -> eltype(t)(0.0)
    regrid_dirpath = joinpath(pkgdir(ClimaLSM), "test", "soil_static")
    raw_datapath = joinpath(pkgdir(ClimaLSM), "means_10.0.nc")
    runoff = ClimaLSM.Soil.TOPMODELRunoff{FT}(
        0.2,
        regrid_dirpath,
        raw_datapath,
        ClimaLSM.Soil.TOPMODELSubsurfaceRunoff{FT}(),
    )
    top_flux_bc = ClimaLSM.Soil.RichardsAtmosDrivenFluxBC(precip, runoff)
    sources = ()
    boundary_fluxes =
        (; top = (water = top_flux_bc,), bottom = (water = bot_flux_bc,))
    params = Soil.RichardsParameters{FT, typeof(hcm)}(ν, hcm, K_sat, S_s, θ_r)

    soil = Soil.RichardsModel{FT}(;
        parameters = params,
        domain = soil_domain,
        boundary_conditions = boundary_fluxes,
        sources = sources,
        lateral_flow = false,
    )
    Y, p, cds = initialize(soil)
    t0 = FT(0)
    set_initial_aux_state! = make_set_initial_aux_state(soil)
    set_initial_aux_state!(p, Y, t0)
    surface_space = soil_domain.space.surface
    comms = surface_space.topology.context
    landsea_mask_manual = regrid_netcdf_to_field(
        FT,
        regrid_dirpath,
        comms,
        raw_datapath,
        "landsea_mask",
        surface_space,
    )

    ϕ̄_manual = regrid_netcdf_to_field(
        FT,
        regrid_dirpath,
        comms,
        input_file,
        "ϕ_mean",
        surface_space,
    )

    fmax_manual = regrid_netcdf_to_field(
        FT,
        regrid_dirpath,
        comms,
        input_file,
        "fmax",
        surface_space,
    )
    @test p.soil.landsea_mask == landsea_mask_manual
    @test p.soil.fmax == fmax_manual
    @test p.soil.ϕ̄ == ϕ̄_manual

end
rm(regrid_dirpath; recursive = true, force = true)

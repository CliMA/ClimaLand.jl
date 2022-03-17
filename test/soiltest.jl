saved_values = SavedValues(FT, ClimaCore.Fields.FieldVector)
ν = FT(0.495);
Ksat = FT(0.0443 / 3600 / 100); # m/s
S_s = FT(1e-3); #inverse meters
vg_n = FT(2.0);
vg_α = FT(2.6); # inverse meters
vg_m = FT(1) - FT(1) / vg_n;
θ_r = FT(0);
zmax = FT(0);
zmin = FT(-10);
nelems = 50;

soil_domain = Column(FT, zlim = (zmin, zmax), nelements = nelems);
top_flux_bc = FT(0.0);
bot_flux_bc = FT(0.0);
sources = ()
boundary_fluxes = FluxBC{FT}(top_flux_bc, bot_flux_bc)
params = Soil.RichardsParameters{FT}(ν, vg_α, vg_n, vg_m, Ksat, S_s, θ_r);

soil = Soil.RichardsModel{FT}(;
    param_set = params,
    domain = soil_domain,
    boundary_conditions = boundary_fluxes,
    sources = sources,
)

Y, p, coords = initialize(soil)

# specify ICs
function init_soil!(Ysoil, z, params)
    function hydrostatic_profile(
        z::FT,
        params::RichardsParameters{FT},
    ) where {FT}
        @unpack ν, vg_α, vg_n, vg_m, θ_r = params
        #unsaturated zone only, assumes water table starts at z_∇
        z_∇ = FT(-10)# matches zmin
        S = FT((FT(1) + (vg_α * (z - z_∇))^vg_n)^(-vg_m))
        ϑ_l = S * (ν - θ_r) + θ_r
        return FT(ϑ_l)
    end
    Ysoil.soil.ϑ_l .= hydrostatic_profile.(z, Ref(params))
end

init_soil!(Y, coords.z, soil.param_set)

soil_ode! = make_ode_function(soil)

t0 = FT(0);
tf = FT(60);
dt = FT(1);
cb = SavingCallback((u, t, integrator) -> copy(integrator.p), saved_values)
prob = ODEProblem(soil_ode!, Y, (t0, tf), p);
sol = solve(prob, Euler(); dt = dt, callback = cb);

@test sum(parent(sol.u[end]) .== parent(Y.soil.ϑ_l)) == nelems
# should be hydrostatic equilibrium at every layer, at each step:
@test mean(
    sum([
        parent(saved_values.saveval[k].soil.ψ .+ coords.z)[:] .+ 10.0 for
        k in 2:1:50
    ]),
) < 1e-10

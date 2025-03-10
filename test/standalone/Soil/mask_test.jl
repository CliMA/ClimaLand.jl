using ClimaComms
apply_mask = true
ClimaComms.@import_required_backends
using ClimaLand
using ClimaLand.Soil
using ClimaLand.Domains
using ClimaComms
using ClimaUtilities.ClimaArtifacts
import Interpolations
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.Regridders: InterpolationsRegridder
using ClimaCore: DataLayouts, Fields, Spaces

FT = Float64
radius = FT(6378.1e3)
depth = FT(50)
domain = ClimaLand.Domains.SphericalShell(;
    radius = radius,
    depth = depth,
    nelements = (101, 15),
    npolynomial = 1,
    dz_tuple = FT.((10.0, 0.05)),
);

# Get the mask and interpolate to our surface space
surface_space = domain.space.surface # 2d space
mask_path = ClimaLand.Artifacts.landseamask_file_path(; resolution = "60arcs")

threshold = 0.5
regridder_type = :InterpolationsRegridder
extrapolation_bc =
    (Interpolations.Periodic(), Interpolations.Flat(), Interpolations.Flat())
context = ClimaComms.context(surface_space)
mask = SpaceVaryingInput(
    mask_path,
    "landsea",
    surface_space;
    regridder_type,
    regridder_kwargs = (; extrapolation_bc,),
)
apply_threshold(field, value) =
    field > value ? eltype(field)(1) : eltype(field)(0)
binary_mask = apply_threshold.(mask, threshold)
if apply_mask
    Spaces.set_mask!(surface_space, binary_mask)
end

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

Y, p, coords = initialize(soil); # Y is the state vector, initially all zero
# test 2d
sfc_field = ClimaLand.top_center_to_surface(Y.soil.ϑ_l);
sfc_field .= 1.0;
@test extrema(sfc_field) == (0.0, 1.0)

# test 3d
Y.soil.ϑ_l .= 1.0
@test extrema(Y.soil.ϑ_l) == (0.0, 1.0)


# More complex function
t0 = FT(0)
set_initial_cache! = make_set_initial_cache(soil)
set_initial_cache!(p, Y, t0)
dY = similar(Y)
@. dY = 0

imp_tendency! = make_imp_tendency(soil)
imp_tendency!(dY, Y, p, t0)

# Test that the masked parts of dY did not update and are still zero
mask_indices = parent(binary_mask) .< eps(FT)
@test extrema(parent(dY.soil.ϑ_l)[:, mask_indices]) == (0.0, 0.0)

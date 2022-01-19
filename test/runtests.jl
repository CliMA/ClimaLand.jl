using Test

using ClimaLSM
using ClimaLSM.Domains: Column, RootDomain
include(joinpath("../src", "Soil.jl"))
using .Soil
include(joinpath("../src", "Roots.jl"))
using .Roots

FT = Float64
soil_domain = Column(FT, zlim = (-1.0, 0.0), nelements = 10)
soil = Soil.RichardsModel{FT}(; param_set = nothing, domain = soil_domain)

root_domain = RootDomain{FT}([-0.5], [0.0, 1.0])
roots = Roots.RootsModel{FT}(; domain = root_domain, param_set = nothing)
# In the future, the user would call initialize and get back coordinates and the state
# structure (empty). They then would initialize Y =Y0 using their IC function(coords)
# then call simulations(model, Y0). P would not need to be shown to the user, as it
# can be created and initialized to the correct starting values (consistent with Y0)
# within Simulations. 
Ysoil, psoil, csoil = initialize(soil)
Yroots, proots, croots = initialize(roots)

# vs
land = LandModel{FT, typeof(soil), typeof(roots)}(soil, roots)
Yland, pland, cland = initialize(land)

# Same structure
@test propertynames(pland.soil) == propertynames(psoil.soil)
@test propertynames(Yland.soil) == propertynames(Ysoil.soil)
@test propertynames(Yland.roots) == propertynames(Yroots.roots)
@test propertynames(pland.roots) == propertynames(proots.roots)

# Make sure this runs...
odef! = make_ode_function(land)
dY = similar(Yland)
odef!(dY, Yland, pland, 0.0)

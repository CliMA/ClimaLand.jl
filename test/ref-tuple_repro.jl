# Initially from soiltest.jl phase change source term test
# To reproduce, run this script on GPU
import ClimaComms
ClimaComms.@import_required_backends
using Statistics
using ClimaCore
import ClimaParams as CP
using ClimaLand
using ClimaLand.Domains: Column
using ClimaLand.Soil
import ClimaLand.Parameters as LP

FT = Float64
earth_param_set = LP.LandParameters(FT)
ν = FT(0.495)
θ_r = FT(0.1)
hydrology_cm = vanGenuchten{FT}(; α = FT(2.6), n = FT(2.0))

soil_domain = Column(; zlim = (FT(-1), FT(0)), nelements = 200)
space_3d = soil_domain.space.subsurface # CenterFiniteDifferenceSpace

θ_l = ClimaCore.Fields.ones(space_3d)
θ_i = ClimaCore.Fields.ones(space_3d)
T = ClimaCore.Fields.ones(space_3d)
κ = ClimaCore.Fields.ones(space_3d)
tau = ClimaCore.Fields.ones(space_3d)

# This fails with dynamic function invocation when `LandParameters`
# and `vanGenuchten` both use `tuple` for broadcasting, and it
# passes when `Ref` is used for either `LandParameters` or `vanGenuchten` broadcasting
@. -phase_change_source(θ_l, θ_i, T, tau, ν, θ_r, hydrology_cm, earth_param_set)


# # Broadcasting f2 fails with the same error as the original matric_potential call
# function f2(hydrology_cm, earth_param_set, ν, θ_r, θ_l)
#     _ρ_i = FT(LP.ρ_cloud_ice(earth_param_set))
#     _ρ_l = FT(LP.ρ_cloud_liq(earth_param_set))
#     _LH_f0 = FT(LP.LH_f0(earth_param_set))
#     _T_freeze = FT(LP.T_freeze(earth_param_set))
#     _grav = FT(LP.grav(earth_param_set))
#     # According to Dall'Amico (text above equation 1), ψw0 corresponds
#     # to the matric potential corresponding to the total water content (liquid and ice).
#     θtot = min(_ρ_i / _ρ_l * θ_i + θ_l, ν)

#     matric_potential(hydrology_cm, effective_saturation(ν, θtot, θ_r))
# end

# @. f2(hydrology_cm, earth_param_set, ν, θ_r, θ_l)


# # Defining θtot either of these ways, then calling f doesn't fail
# function f(hydrology_cm, ν, θtot, θ_r)
#     @. matric_potential(hydrology_cm, effective_saturation(ν, θtot, θ_r))
# end

# _ρ_i = FT(LP.ρ_cloud_ice(earth_param_set))
# _ρ_l = FT(LP.ρ_cloud_liq(earth_param_set))
# _LH_f0 = FT(LP.LH_f0(earth_param_set))
# _T_freeze = FT(LP.T_freeze(earth_param_set))
# _grav = FT(LP.grav(earth_param_set))
# θtot = @. min(_ρ_i / _ρ_l * θ_i + θ_l, ν) # f returns [-0.0, -0.0, ..., -0.0]

# # θtot = ClimaCore.Fields.ones(space_3d) # f returns [NaN, NaN, ..., NaN]
# f(hydrology_cm, ν, θtot, θ_r)

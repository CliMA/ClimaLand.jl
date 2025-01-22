module Parameters

abstract type AbstractLandParameters end
const ALP = AbstractLandParameters

Base.@kwdef struct LandParameters{FT, TP, SFP, IP} <: ALP
    K_therm::FT
    ρ_cloud_liq::FT
    ρ_cloud_ice::FT
    cp_l::FT
    cp_i::FT
    T_0::FT
    LH_v0::FT
    LH_s0::FT
    Stefan::FT
    T_freeze::FT
    grav::FT
    MSLP::FT
    D_vapor::FT
    gas_constant::FT
    molmass_water::FT
    h_Planck::FT
    light_speed::FT
    avogad::FT
    thermo_params::TP
    surf_flux_params::SFP
    insol_params::IP
end

Base.eltype(::LandParameters{FT}) where {FT} = FT
Base.broadcastable(ps::LandParameters) = tuple(ps)

# wrapper methods:
P_ref(ps::ALP) = ps.MSLP
K_therm(ps::ALP) = ps.K_therm
ρ_cloud_liq(ps::ALP) = ps.ρ_cloud_liq
ρ_cloud_ice(ps::ALP) = ps.ρ_cloud_ice
cp_l(ps::ALP) = ps.cp_l
cp_i(ps::ALP) = ps.cp_i
T_0(ps::ALP) = ps.T_0
LH_v0(ps::ALP) = ps.LH_v0
LH_s0(ps::ALP) = ps.LH_s0
Stefan(ps::ALP) = ps.Stefan
T_freeze(ps::ALP) = ps.T_freeze
grav(ps::ALP) = ps.grav
D_vapor(ps::ALP) = ps.D_vapor
gas_constant(ps::ALP) = ps.gas_constant
molar_mass_water(ps::ALP) = ps.molmass_water
planck_constant(ps::ALP) = ps.h_Planck
avogadro_constant(ps::ALP) = ps.avogad
light_speed(ps::ALP) = ps.light_speed
# Derived parameters
LH_f0(ps::ALP) = LH_s0(ps) - LH_v0(ps)
ρ_m_liq(ps::ALP) = ρ_cloud_liq(ps) / molar_mass_water(ps)
# Dependency parameter wrappers
thermodynamic_parameters(ps::ALP) = ps.thermo_params
surface_fluxes_parameters(ps::ALP) = ps.surf_flux_params
insolation_parameters(ps::ALP) = ps.insol_params

end # module

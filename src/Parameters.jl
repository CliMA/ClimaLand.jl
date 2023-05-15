module Parameters

abstract type AbstractLSMParameters end
const ALSMP = AbstractLSMParameters

Base.@kwdef struct LSMParameters{FT, TP, SFP, IP} <: ALSMP
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

# wrapper methods:
P_ref(ps::ALSMP) = ps.MSLP
K_therm(ps::ALSMP) = ps.K_therm
ρ_cloud_liq(ps::ALSMP) = ps.ρ_cloud_liq
ρ_cloud_ice(ps::ALSMP) = ps.ρ_cloud_ice
cp_l(ps::ALSMP) = ps.cp_l
cp_i(ps::ALSMP) = ps.cp_i
T_0(ps::ALSMP) = ps.T_0
LH_v0(ps::ALSMP) = ps.LH_v0
LH_s0(ps::ALSMP) = ps.LH_s0
Stefan(ps::ALSMP) = ps.Stefan
T_freeze(ps::ALSMP) = ps.T_freeze
grav(ps::ALSMP) = ps.grav
D_vapor(ps::ALSMP) = ps.D_vapor
gas_constant(ps::ALSMP) = ps.gas_constant
molar_mass_water(ps::ALSMP) = ps.molmass_water
planck_constant(ps::ALSMP) = ps.h_Planck
avogadro_constant(ps::ALSMP) = ps.avogad
light_speed(ps::ALSMP) = ps.light_speed
# Derived parameters
LH_f0(ps::ALSMP) = LH_s0(ps) - LH_v0(ps)
ρ_m_liq(ps::ALSMP) = ρ_cloud_liq(ps) / molar_mass_water(ps)
# Dependency parameter wrappers
thermodynamic_parameters(ps::ALSMP) = ps.thermo_params
surface_fluxes_parameters(ps::ALSMP) = ps.surf_flux_params
insolation_parameters(ps::ALSMP) = ps.insol_params

end # module

export make_parameters, SiteParams, HeteroRespP, PlantHydraulicsP

"""
    make_parameters(;
        params,
        )

Setup parameters for a site
"""
function make_parameters(site_ID)
    default_params_dict = Dict(
        "US-MOz" => ozark_default_params(),
        "US-Ha1" => harvard_default_params(), # need to add those once I create new harvard_default_params() etc
        "US-NR1" => niwotridge_default_params(),
        "US-Var" => vairaranch_default_params(),
    )

    if haskey(default_params_dict, site_ID)
        return default_params_dict[site_ID]
    else
        println("This site_ID does not exist")
    end
end

# The structs below are necessary and different from ClimaLand/src struct 
# because of, for example, shared parameters (e.g., porosity in soil or soilco2)
# or, another example, parameters calculated from other parameters (e.g., κ_solid)
# the structs below only contain unique params, that are given as number (not computed FTwhere)

struct HeteroRespP # differs from SoilCO2ModelParameters because of porosity
    θ_a100::FT
    D_ref::FT
    b::FT
    D_liq::FT
    α_sx::FT
    Ea_sx::FT
    kM_sx::FT
    kM_o2::FT
    O2_a::FT
    D_oa::FT
    p_sx::FT
end

struct PlantHydraulicsP # PlantHydraulicsParameters has porosity, storativity, which are defined in soil, and conductivity_model, retention_model, which need to be computed
    SAI::FT
    f_root_to_shoot::FT
    capacity::FT
    S_s::FT
    rooting_depth::FT
    conductivity_model::Any
    retention_model::Any
end

struct SiteParams
    hetero_resp::HeteroRespP
    auto_resp::AutotrophicRespirationParameters
    soil::EnergyHydrologyParameters
    radiative_transfer::TwoStreamParameters
    canopy_energy_balance::BigLeafEnergyParameters
    conductance::MedlynConductanceParameters
    photosynthesis::FarquharParameters
    plant_hydraulics::PlantHydraulicsP
end

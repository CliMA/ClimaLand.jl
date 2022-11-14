import CLIMAParameters as CP
import Thermodynamics.Parameters as TDP
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions as UF
import ClimaLSM.Parameters as LSMP

#=
import ClimaLSM
include(joinpath(pkgdir(ClimaLSM), "parameters", "create_parameters.jl"))
param_set = create_lsm_parameters(FT)
=#
function create_lsm_parameters(FT)
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")

    aliases = string.(fieldnames(TDP.ThermodynamicsParameters))
    param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
    thermo_params = TDP.ThermodynamicsParameters{FT}(; param_pairs...)
    TP = typeof(thermo_params)

    aliases = [
        "Pr_0_Businger",
        "a_m_Businger",
        "a_h_Businger",
        "ζ_a_Businger",
        "γ_Businger",
    ]
    pairs = CP.get_parameter_values!(toml_dict, aliases, "UniversalFunctions")
    pairs = (; pairs...) # convert to NamedTuple
    pairs = (;
        Pr_0 = pairs.Pr_0_Businger,
        a_m = pairs.a_m_Businger,
        a_h = pairs.a_h_Businger,
        ζ_a = pairs.ζ_a_Businger,
        γ = pairs.γ_Businger,
    )
    ufp = UF.BusingerParams{FT}(; pairs...)
    UFP = typeof(ufp)

    pairs = CP.get_parameter_values!(
        toml_dict,
        ["von_karman_const"],
        "SurfaceFluxesParameters",
    )
    surf_flux_params =
        SFP.SurfaceFluxesParameters{FT, UFP, TP}(; pairs..., ufp, thermo_params)
    SFPS = typeof(surf_flux_params)

    aliases = string.(fieldnames(LSMP.LSMParameters))
    param_pairs = CP.get_parameter_values!(toml_dict, aliases, "ClimaLSM")
    param_set = LSMP.LSMParameters{FT, TP, SFPS}(;
        param_pairs...,
        thermo_params,
        surf_flux_params,
    )

    # logfilepath = joinpath(@__DIR__, "logfilepath_$FT.toml")
    # CP.log_parameter_information(toml_dict, logfilepath)
    return param_set
end

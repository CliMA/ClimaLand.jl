if !("." in LOAD_PATH) # for ease of local testing
    push!(LOAD_PATH, ".")
end

import ClimaParams as CP
import SurfaceFluxes
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions
using Thermodynamics
import SurfaceFluxes.Parameters.SurfaceFluxesParameters
import Plots

FT = Float32
param_set = SurfaceFluxesParameters(FT, UniversalFunctions.BusingerParams)
thermo_params = param_set.thermo_params
uft = UniversalFunctions.universal_func_type(typeof(param_set.ufp))

# Define surface parameters. Note that recovery parameters purely depend on LMO, scale variables and Δ(interior - surface) 
# for each variable
θ_sfc = FT(295)
p_sfc = FT(100000)
p_in = FT(99100)

bonan_data_dθz = [
    -3.552713678800501e-15 19.193548387096776
    -0.4848484848484844 19.354838709677423
    -1.90909090909091 19.67741935483872
    -2.6363636363636402 20.48387096774194
    -3.242424242424242 22.258064516129036
    -3.575757575757578 25.161290322580648
    -3.6969696969696972 29.193548387096776
    -3.8787878787878824 38.387096774193544
    -3.9696969696969724 43.87096774193547
    -4 49.35483870967741
]
bonan_data_duz = [
    0 19.67741935483872
    0.39393939393939414 20.000000000000004
    0.8636363636363636 20.967741935483875
    1.2575757575757573 22.419354838709683
    1.5909090909090908 24.67741935483871
    1.8333333333333335 27.741935483870968
    2.0303030303030303 32.41935483870968
    2.1666666666666665 37.096774193548384
    2.287878787878788 42.419354838709666
    2.409090909090909 49.51612903225805
]

data_Δθ = bonan_data_dθz[:, 1]
data_Δu = bonan_data_duz[:, 1]
data_uz = bonan_data_duz[:, 2]
data_θz = bonan_data_dθz[:, 2]

il = size(bonan_data_duz)[2]

LMO = FT(-10)
ustar = FT(0.4)
θstar = FT(-0.5)
canopy_disp = FT(19)
result_u = [];
result_θ = [];

for il in 1:10
    ts_int_test = Thermodynamics.PhaseDry_pθ(
        thermo_params,
        p_in,
        FT(θ_sfc + bonan_data_dθz[il, 1]),
    )
    ts_sfc_test = Thermodynamics.PhaseDry_pθ(thermo_params, p_sfc, θ_sfc)
    sc = SF.ValuesOnly(
        SF.StateValues(
            FT(bonan_data_duz[il, 2] - canopy_disp),
            (FT(bonan_data_duz[il, 1]), FT(0)),
            ts_int_test,
        ),
        SF.StateValues(FT(0), (FT(0), FT(0)), ts_sfc_test),
        FT(0.6),
        FT(0.0816),
    )
    push!(
        result_u,
        SurfaceFluxes.recover_profile(
            param_set,
            sc,
            LMO,
            data_uz[il] - canopy_disp,
            ustar,
            FT(0),
            UniversalFunctions.MomentumTransport(),
            SurfaceFluxes.PointValueScheme(),
        ),
    )
    push!(
        result_θ,
        SurfaceFluxes.recover_profile(
            param_set,
            sc,
            LMO,
            data_uz[il] - canopy_disp,
            θstar,
            θ_sfc,
            UniversalFunctions.HeatTransport(),
            SurfaceFluxes.PointValueScheme(),
        ),
    )
end
pl = Plots.@layout [a b; c d]
p1 = Plots.plot(
    result_u,
    data_uz,
    xlim = (0, 4),
    ylim = (15, 50),
    ls = :dashdot,
    label = "LMO=-10",
)
p2 = Plots.plot(
    result_θ .- θ_sfc,
    data_uz,
    xlim = (-8, 0),
    ylim = (15, 50),
    ls = :dashdot,
    label = "LMO=-10",
)

LMO = FT(-50)
result_u = [];
result_θ = [];
for il in 1:10
    ts_int_test = Thermodynamics.PhaseDry_pθ(
        thermo_params,
        p_in,
        FT(θ_sfc + bonan_data_dθz[il, 1]),
    )
    ts_sfc_test = Thermodynamics.PhaseDry_pθ(thermo_params, p_sfc, θ_sfc)
    sc = SF.ValuesOnly(
        SF.StateValues(
            FT(bonan_data_duz[il, 2] - canopy_disp),
            (FT(bonan_data_duz[il, 1]), FT(0)),
            ts_int_test,
        ),
        SF.StateValues(FT(0), (FT(0), FT(0)), ts_sfc_test),
        FT(0.6),
        FT(0.0816),
    )
    push!(
        result_u,
        SurfaceFluxes.recover_profile(
            param_set,
            sc,
            LMO,
            data_uz[il] - canopy_disp,
            ustar,
            FT(0),
            UniversalFunctions.MomentumTransport(),
            SurfaceFluxes.PointValueScheme(),
        ),
    )
    push!(
        result_θ,
        SurfaceFluxes.recover_profile(
            param_set,
            sc,
            LMO,
            data_uz[il] - canopy_disp,
            θstar,
            θ_sfc,
            UniversalFunctions.HeatTransport(),
            SurfaceFluxes.PointValueScheme(),
        ),
    )
end
p1 = Plots.plot!(
    p1,
    result_u,
    data_uz,
    xlim = (0, 4),
    ylim = (15, 50),
    ls = :dash,
    label = "LMO=-50",
)
p2 = Plots.plot!(
    p2,
    result_θ .- θ_sfc,
    data_uz,
    xlim = (-8, 0),
    ylim = (15, 50),
    ls = :dash,
    label = "LMO=-50",
)

LMO = FT(-1000)
result_u = [];
result_θ = [];
for il in 1:10
    ts_int_test = Thermodynamics.PhaseDry_pθ(
        thermo_params,
        p_in,
        FT(θ_sfc + bonan_data_dθz[il, 1]),
    )
    ts_sfc_test = Thermodynamics.PhaseDry_pθ(thermo_params, p_sfc, θ_sfc)
    sc = SF.ValuesOnly(
        SF.StateValues(
            FT(bonan_data_duz[il, 2] - canopy_disp),
            (FT(bonan_data_duz[il, 1]), FT(0)),
            ts_int_test,
        ),
        SF.StateValues(FT(0), (FT(0), FT(0)), ts_sfc_test),
        FT(0.6),
        FT(0.0816),
    )
    push!(
        result_u,
        SurfaceFluxes.recover_profile(
            param_set,
            sc,
            LMO,
            data_uz[il] - canopy_disp,
            ustar,
            FT(0),
            UniversalFunctions.MomentumTransport(),
            SurfaceFluxes.PointValueScheme(),
        ),
    )
    push!(
        result_θ,
        SurfaceFluxes.recover_profile(
            param_set,
            sc,
            LMO,
            data_uz[il] - canopy_disp,
            θstar,
            θ_sfc,
            UniversalFunctions.HeatTransport(),
            SurfaceFluxes.PointValueScheme(),
        ),
    )
end
p1 = Plots.plot!(
    p1,
    result_u,
    data_uz,
    xlim = (0, 4),
    ylim = (15, 50),
    ls = :solid,
    label = "LMO=-1000",
    xlabel = "u(z) [ms⁻¹]",
    ylabel = "Height, z [m]",
)
p2 = Plots.plot!(
    p2,
    result_θ .- θ_sfc,
    data_uz,
    xlim = (-8, 0),
    ylim = (15, 50),
    ls = :solid,
    label = "LMO=-1000",
    xlabel = "θ(z)-θ_sfc [K]",
    ylabel = "Height, z [m]",
)

LMO = FT(30)
ustar = FT(0.13)
θstar = FT(0.06)
canopy_disp = FT(19)
result_u = [];
result_θ = [];

for il in 1:10
    ts_int_test = Thermodynamics.PhaseDry_pθ(
        thermo_params,
        p_in,
        FT(θ_sfc + bonan_data_dθz[il, 1]),
    )
    ts_sfc_test = Thermodynamics.PhaseDry_pθ(thermo_params, p_sfc, θ_sfc)
    sc = SF.ValuesOnly(
        SF.StateValues(
            FT(bonan_data_duz[il, 2] - canopy_disp),
            (FT(bonan_data_duz[il, 1]), FT(0)),
            ts_int_test,
        ),
        SF.StateValues(FT(0), (FT(0), FT(0)), ts_sfc_test),
        FT(0.6),
        FT(0.0816),
    )
    push!(
        result_u,
        SurfaceFluxes.recover_profile(
            param_set,
            sc,
            LMO,
            data_uz[il] - canopy_disp,
            ustar,
            FT(0),
            UniversalFunctions.MomentumTransport(),
            SurfaceFluxes.PointValueScheme(),
        ),
    )
    push!(
        result_θ,
        SurfaceFluxes.recover_profile(
            param_set,
            sc,
            LMO,
            data_uz[il] - canopy_disp,
            θstar,
            θ_sfc,
            UniversalFunctions.HeatTransport(),
            SurfaceFluxes.PointValueScheme(),
        ),
    )
end
p3 = Plots.plot(
    result_u,
    data_uz,
    xlim = (0, 4),
    ylim = (15, 50),
    ls = :dashdot,
    label = "LMO=30",
)
p4 = Plots.plot(
    result_θ .- θ_sfc,
    data_uz,
    xlim = (0, 2),
    ylim = (15, 50),
    ls = :dashdot,
    label = "LMO=30",
)

LMO = FT(50)
result_u = [];
result_θ = [];
for il in 1:10
    ts_int_test = Thermodynamics.PhaseDry_pθ(
        thermo_params,
        p_in,
        FT(θ_sfc + bonan_data_dθz[il, 1]),
    )
    ts_sfc_test = Thermodynamics.PhaseDry_pθ(thermo_params, p_sfc, θ_sfc)
    sc = SF.ValuesOnly(
        SF.StateValues(
            FT(bonan_data_duz[il, 2] - canopy_disp),
            (FT(bonan_data_duz[il, 1]), FT(0)),
            ts_int_test,
        ),
        SF.StateValues(FT(0), (FT(0), FT(0)), ts_sfc_test),
        FT(0.6),
        FT(0.0816),
    )
    push!(
        result_u,
        SurfaceFluxes.recover_profile(
            param_set,
            sc,
            LMO,
            data_uz[il] - canopy_disp,
            ustar,
            FT(0),
            UniversalFunctions.MomentumTransport(),
            SurfaceFluxes.PointValueScheme(),
        ),
    )
    push!(
        result_θ,
        SurfaceFluxes.recover_profile(
            param_set,
            sc,
            LMO,
            data_uz[il] - canopy_disp,
            θstar,
            θ_sfc,
            UniversalFunctions.HeatTransport(),
            SurfaceFluxes.PointValueScheme(),
        ),
    )
end
p3 = Plots.plot!(
    p3,
    result_u,
    data_uz,
    xlim = (0, 4),
    ylim = (15, 50),
    ls = :dash,
    label = "LMO=50",
)
p4 = Plots.plot!(
    p4,
    result_θ .- θ_sfc,
    data_uz,
    xlim = (0, 2),
    ylim = (15, 50),
    ls = :dash,
    label = "LMO=50",
)

LMO = FT(1000)
result_u = [];
result_θ = [];
for il in 1:10
    ts_int_test = Thermodynamics.PhaseDry_pθ(
        thermo_params,
        p_in,
        FT(θ_sfc + bonan_data_dθz[il, 1]),
    )
    ts_sfc_test = Thermodynamics.PhaseDry_pθ(thermo_params, p_sfc, θ_sfc)
    sc = SF.ValuesOnly(
        SF.StateValues(
            FT(bonan_data_duz[il, 2] - canopy_disp),
            (FT(bonan_data_duz[il, 1]), FT(0)),
            ts_int_test,
        ),
        SF.StateValues(FT(0), (FT(0), FT(0)), ts_sfc_test),
        FT(0.6),
        FT(0.0816),
    )
    push!(
        result_u,
        SurfaceFluxes.recover_profile(
            param_set,
            sc,
            LMO,
            data_uz[il] - canopy_disp,
            ustar,
            FT(0),
            UniversalFunctions.MomentumTransport(),
            SurfaceFluxes.PointValueScheme(),
        ),
    )
    push!(
        result_θ,
        SurfaceFluxes.recover_profile(
            param_set,
            sc,
            LMO,
            data_uz[il] - canopy_disp,
            θstar,
            θ_sfc,
            UniversalFunctions.HeatTransport(),
            SurfaceFluxes.PointValueScheme(),
        ),
    )
end
p3 = Plots.plot!(
    p3,
    result_u,
    data_uz,
    xlim = (0, 4),
    ylim = (15, 50),
    ls = :solid,
    label = "LMO=1000",
    xlabel = "u(z) [ms⁻¹]",
    ylabel = "Height, z [m]",
)
p4 = Plots.plot!(
    p4,
    result_θ .- θ_sfc,
    data_uz,
    xlim = (0, 2),
    ylim = (15, 50),
    ls = :solid,
    label = "LMO=1000",
    xlabel = "θ(z)-θ_sfc [K]",
    ylabel = "Height, z [m]",
)

bonan_fig = Plots.plot(p1, p2, p3, p4; layout = pl)
Plots.savefig("Bonan_Fig6-4.svg")

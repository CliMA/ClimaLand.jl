import NCDatasets as NC
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF
using Statistics
import Thermodynamics
import ArtifactWrappers
const AW = ArtifactWrappers
const TD = Thermodynamics

FT = Float32
param_set = SFP.SurfaceFluxesParameters(FT, BusingerParams)
thermo_params = param_set.thermo_params
uft = UF.BusingerType()

#! format: off
PyCLES_output_dataset = AW.ArtifactWrapper(
    @__DIR__,
    "PyCLES_output",
    AW.ArtifactFile[
    AW.ArtifactFile(url = "https://caltech.box.com/shared/static/johlutwhohvr66wn38cdo7a6rluvz708.nc", filename = "Rico.nc",),
    AW.ArtifactFile(url = "https://caltech.box.com/shared/static/zraeiftuzlgmykzhppqwrym2upqsiwyb.nc", filename = "Gabls.nc",),
    AW.ArtifactFile(url = "https://caltech.box.com/shared/static/toyvhbwmow3nz5bfa145m5fmcb2qbfuz.nc", filename = "DYCOMS_RF01.nc",),
    AW.ArtifactFile(url = "https://caltech.box.com/shared/static/jci8l11qetlioab4cxf5myr1r492prk6.nc", filename = "Bomex.nc",),
    ],
)
#! format: on
PyCLES_output_dataset_path = AW.get_data_folder(PyCLES_output_dataset)

files = ["DYCOMS_RF01.nc", "Bomex.nc", "Rico.nc", "Gabls.nc"];
for f in files
    @info "Casename: $f"

    nt = NC.NCDataset(joinpath(PyCLES_output_dataset_path, f), "r") do data
        prof = data.group["profiles"]
        timeseries = data.group["timeseries"]
        z = Array(data["z_half"])
        u = Array(prof["u_mean"])
        v = Array(prof["v_mean"])
        qt = Array(prof["qt_mean"])
        ql = Array(prof["ql_mean"])
        ρ = Array(prof["rho"])
        p0 = Array(prof["p0"])
        b = Array(prof["buoyancy_mean"])
        θli = Array(prof["thetali_mean"]) # Care with variable names across cases ()?
        T = Array(prof["temperature_mean"]) # Care with variable names across cases ()?
        SHF = Array(timeseries["shf_surface_mean"]) # Care with variable names across cases ()?
        LHF = Array(timeseries["lhf_surface_mean"]) # Care with variable names across cases ()?
        (; z, u, v, qt, ql, ρ, p0, b, θli, T, SHF, LHF)
    end
    z, u, v, qt, ql, ρ, p0, b, θli, T, SHF, LHF = nt

    function getval(X) # Consider only the last timestep
        X[:, end]
    end

    # Data at first interior node (x_in)
    # TODO Make sure that the first node ii = 1 is at the "surface"
    # for the tests to be consistent
    ii = 2

    shf = mean(getval(SHF))
    lhf = mean(getval(LHF))

    z_in = z[ii]
    u_in = getval(u)[ii]
    v_in = getval(v)[ii]
    b_in = getval(b)[ii]
    ρ_in = getval(ρ)[ii]
    θ_in = getval(θli)[ii]
    T_in = getval(T)[ii]
    qt_in = getval(qt)[ii]

    # Surface values for variables
    jj = 1
    u_sfc = FT(0)
    v_sfc = FT(0)
    b_sfc = getval(b)[jj]
    ρ_sfc = getval(ρ)[jj]
    qt_sfc = getval(qt)[jj]
    θ_sfc = getval(θli)[jj]
    T_sfc = getval(T)[jj]
    z_sfc = FT(0)

    # Roughness
    ## Roughness lengths
    z0m = FT(0.001)
    z0b = FT(0.001)

    ts_sfc = TD.PhaseEquil_ρθq(thermo_params, ρ_sfc, θ_sfc, qt_sfc)
    ts_in = TD.PhaseEquil_ρθq(thermo_params, ρ_in, θ_in, qt_in)

    u_in = (FT(u_in), FT(v_in))
    u_sfc = (FT(u_sfc), FT(v_sfc))

    state_sfc = SF.StateValues(z_sfc, u_sfc, ts_sfc)
    state_in = SF.StateValues(z_in, u_in, ts_in)

    kwargs = (; state_in, state_sfc, z0m, z0b)

    # shf and lhf in dycoms are mising /wrong in the file
    if f == "DYCOMS_RF01.nc"
        sc = SF.Fluxes(state_in, state_sfc, FT(15), FT(115), z0m, z0b)
    elseif f == "Bomex.nc"
        sc = SF.FluxesAndFrictionVelocity(
            state_in,
            state_sfc,
            shf,
            lhf,
            FT(0.28),
            z0m,
            z0b,
        )
    elseif f == "Rico.nc"
        sc = SF.Coefficients(
            state_in,
            state_sfc,
            FT(0.001229),
            FT(0.001094),
            z0m,
            z0b,
        )
    elseif f == "Gabls.nc"
        sc = SF.ValuesOnly(state_in, state_sfc, z0m, z0b)
    end
    result = SF.surface_conditions(param_set, sc)
end

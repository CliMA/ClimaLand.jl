#!/bin/bash

set -euo pipefail

cd "$(dirname "$0")"

PASS_COUNT=0
FAIL_COUNT=0
JULIA_PROJECT=".buildkite"
JULIA_BIN="${JULIA_BIN:-/resnick/groups/esm/software/julia/julia-1.12.5/bin/julia}"

pass() {
    echo "[PASS] $1"
    PASS_COUNT=$((PASS_COUNT + 1))
}

fail() {
    echo "[FAIL] $1"
    FAIL_COUNT=$((FAIL_COUNT + 1))
}

check_file() {
    local path="$1"
    if [[ -f "$path" ]]; then
        pass "file exists: $path"
    else
        fail "missing file: $path"
    fi
}

check_bash_syntax() {
    local path="$1"
    if bash -n "$path"; then
        pass "bash syntax: $path"
    else
        fail "bash syntax error: $path"
    fi
}

check_grep() {
    local pattern="$1"
    local path="$2"
    local description="$3"
    if grep -Eq "$pattern" "$path"; then
        pass "$description"
    else
        fail "$description"
    fi
}

check_not_grep() {
    local pattern="$1"
    local path="$2"
    local description="$3"
    if grep -Eq "$pattern" "$path"; then
        fail "$description"
    else
        pass "$description"
    fi
}

check_julia_smoke() {
    local description="$1"
    local script="$2"
    if [[ ! -x "$JULIA_BIN" ]]; then
        fail "$description"
        echo "Julia binary not found: $JULIA_BIN"
        return
    fi
    if "$JULIA_BIN" --project="$JULIA_PROJECT" --startup-file=no --history-file=no -e "$script" >/tmp/callmip_julia_smoke.log 2>&1; then
        pass "$description"
    else
        fail "$description"
        echo "See /tmp/callmip_julia_smoke.log"
    fi
}

echo "=== CalLMIP Pipeline Readiness Test ==="
echo "Repo root: $PWD"
echo

echo "[1/6] Checking critical files"
check_file "submit_callmip_pipeline.sh"
check_file "experiments/calibrate_dk_sor/slurm_calibration.sh"
check_file "experiments/callmip_uq_dk_sor/slurm_emulate_sample.sh"
check_file "experiments/callmip_uq_dk_sor/slurm_callmip_simulations.sh"
check_file "experiments/callmip_uq_dk_sor/slurm_postprocess.sh"
check_file "experiments/calibrate_dk_sor/generate_observations.jl"
check_file "experiments/calibrate_dk_sor/run_calibration.jl"
check_file "experiments/callmip_uq_dk_sor/emulate_sample.jl"
check_file "experiments/callmip_uq_dk_sor/run_callmip_simulations.jl"
check_file "experiments/callmip_uq_dk_sor/write_callmip_netcdf.jl"
check_file "experiments/calibrate_dk_sor/model_interface.jl"
check_file "experiments/calibrate_dk_sor/run_alexis_test.jl"
check_file "experiments/callmip_uq_dk_sor/run_posterior_ensemble.jl"
check_file "ext/fluxnet_simulations/DK-Sor.jl"
check_file "DK_Sor/DK-Sor_daily_aggregated_1997-2013_FLUXNET2015_Flux.nc"
check_file "DK_Sor/DK-Sor_1997-2014_FLUXNET2015_Met.nc"

echo
echo "[2/6] Checking script syntax"
check_bash_syntax "submit_callmip_pipeline.sh"
check_bash_syntax "experiments/calibrate_dk_sor/slurm_calibration.sh"
check_bash_syntax "experiments/callmip_uq_dk_sor/slurm_emulate_sample.sh"
check_bash_syntax "experiments/callmip_uq_dk_sor/slurm_callmip_simulations.sh"
check_bash_syntax "experiments/callmip_uq_dk_sor/slurm_postprocess.sh"

echo
echo "[3/6] Checking observation artifact schema compatibility"
check_grep "y_obs =" "experiments/calibrate_dk_sor/generate_observations.jl" "generate_observations writes calibration y_obs key"
check_grep "noise_cov =" "experiments/calibrate_dk_sor/generate_observations.jl" "generate_observations writes CES noise_cov key"
check_grep "obs_dates =" "experiments/calibrate_dk_sor/generate_observations.jl" "generate_observations writes CES obs_dates key"
check_grep "window_names =|window_dates =|window_pairs =" "experiments/calibrate_dk_sor/generate_observations.jl" "generate_observations writes minibatch window metadata"

if [[ -f "experiments/calibrate_dk_sor/observations.jld2" ]]; then
    echo "[WARN] observations.jld2 exists; stage 1 will regenerate it before calibration"
else
    echo "[WARN] observations.jld2 not present yet; stage 1 will generate it"
fi

echo
echo "[4/6] Checking stage wiring expectations"
check_grep "analyze_posterior_ensemble.jl" "experiments/callmip_uq_dk_sor/slurm_postprocess.sh" "postprocess includes posterior analysis"
check_grep "evaluate_calibration.jl" "experiments/callmip_uq_dk_sor/slurm_postprocess.sh" "postprocess includes calibration evaluation"
check_grep "write_callmip_netcdf.jl" "experiments/callmip_uq_dk_sor/slurm_postprocess.sh" "postprocess includes NetCDF writer"
check_grep "DateTime\(1996, 1, 1\)" "experiments/callmip_uq_dk_sor/callmip_model_interface.jl" "callmip forward model includes 1996 spinup start"
check_grep "OUTPUT_START_DATE = Date\(1997, 1, 1\)" "experiments/callmip_uq_dk_sor/callmip_model_interface.jl" "callmip output starts in 1997"
check_grep "const DT[[:space:]]*= Float64\(900\)" "experiments/callmip_uq_dk_sor/run_callmip_simulations.jl" "callmip simulation DT matches calibration"
check_grep "const DT[[:space:]]*= Float64\(900\)" "experiments/callmip_uq_dk_sor/run_posterior_ensemble.jl" "posterior ensemble DT matches calibration"
check_grep "DKSorModelInterface\(\)" "experiments/calibrate_dk_sor/run_calibration.jl" "calibration instantiates the DK-Sor interface"
check_grep "ClimaCalibrate\.WorkerBackend\(\)" "experiments/calibrate_dk_sor/run_calibration.jl" "calibration uses WorkerBackend()"
check_grep "model_interface," "experiments/calibrate_dk_sor/run_calibration.jl" "calibration passes the interface object to ClimaCalibrate.calibrate"
check_grep "include\(\"fluxnet_simulations/DK-Sor.jl\"\)" "ext/FluxnetSimulationsExt.jl" "Fluxnet extension registers the DK-Sor site file"
check_grep "function FluxnetSimulations\.prescribed_forcing_netcdf" "ext/fluxnet_simulations/forcing.jl" "Fluxnet extension defines prescribed_forcing_netcdf"
check_grep "ClimaCalibrate\.forward_model\(model_interface, 999, 1\)" "experiments/calibrate_dk_sor/run_alexis_test.jl" "Alexis smoke test uses the DK-Sor interface object"
check_grep "ClimaCalibrate\.forward_model\(model_interface, 0, m\)" "experiments/callmip_uq_dk_sor/run_posterior_ensemble.jl" "posterior ensemble uses the DK-Sor interface object"
check_not_grep "n_stem[[:space:]]*=|h_stem[[:space:]]*=|h_leaf[[:space:]]*=" "experiments/calibrate_dk_sor/model_interface.jl" "calibration interface avoids deprecated PlantHydraulicsModel kwargs"
check_not_grep "n_stem[[:space:]]*=|h_stem[[:space:]]*=|h_leaf[[:space:]]*=" "experiments/callmip_uq_dk_sor/callmip_model_interface.jl" "callmip interface avoids deprecated PlantHydraulicsModel kwargs"
check_not_grep "n_stem[[:space:]]*=|h_stem[[:space:]]*=|h_leaf[[:space:]]*=" "experiments/calibrate_dk_sor/run_prior_mean.jl" "prior-mean script avoids deprecated PlantHydraulicsModel kwargs"
check_not_grep "hydraulics\.n_stem|hydraulics\.n_leaf|hydraulics\.ϑ_l\.:\(\$i\)" "experiments/calibrate_dk_sor/model_interface.jl" "calibration interface uses modern hydraulics state initialization"
check_not_grep "hydraulics\.n_stem|hydraulics\.n_leaf|hydraulics\.ϑ_l\.:\(\$i\)" "experiments/callmip_uq_dk_sor/callmip_model_interface.jl" "callmip interface uses modern hydraulics state initialization"
check_not_grep "hydraulics\.n_stem|hydraulics\.n_leaf|hydraulics\.ϑ_l\.:\(\$i\)" "experiments/calibrate_dk_sor/run_prior_mean.jl" "prior-mean script uses modern hydraulics state initialization"
check_not_grep "canopy\.energy\.T" "src/standalone/Vegetation/autotrophic_respiration.jl" "autotrophic respiration avoids stale canopy energy temperature field"
check_grep "p\.drivers\.T" "src/standalone/Vegetation/autotrophic_respiration.jl" "autotrophic respiration uses driver air temperature"

if [[ -x "$JULIA_BIN" ]]; then
    check_julia_smoke \
        "Julia smoke test: DK-Sor Fluxnet site methods exist" \
        'import ClimaLand; import ClimaLand.FluxnetSimulations as FS; FT = Float64; site = Val(:DK_Sor); @assert hasmethod(FS.get_domain_info, Tuple{DataType, typeof(site)}); @assert hasmethod(FS.get_location, Tuple{DataType, typeof(site)}); @assert hasmethod(FS.get_fluxtower_height, Tuple{DataType, typeof(site)}); d = FS.get_domain_info(FT, site); loc = FS.get_location(FT, site); tower = FS.get_fluxtower_height(FT, site); @assert d.nelements == 20; @assert d.zmin < 0; @assert loc.time_offset == 1; @assert tower.atmos_h > 0; println("DK-Sor site-method smoke test passed")'
    check_julia_smoke \
        "Julia smoke test: DK-Sor interface loads and exposes required methods" \
        'include("experiments/calibrate_dk_sor/model_interface.jl"); using Main: DKSorModelInterface; @assert hasmethod(ClimaCalibrate.model_interface_filepath, Tuple{DKSorModelInterface}); @assert hasmethod(ClimaCalibrate.experiment_dir, Tuple{DKSorModelInterface}); @assert hasmethod(ClimaCalibrate.forward_model, Tuple{DKSorModelInterface, Int, Int}); @assert hasmethod(ClimaCalibrate.observation_map, Tuple{DKSorModelInterface, Int}); println("DK-Sor interface smoke test passed")'
    check_julia_smoke \
        "Julia smoke test: NetCDF forcing symbol is available" \
        'import ClimaLand; import ClimaLand.FluxnetSimulations as FS; @assert isdefined(FS, :prescribed_forcing_netcdf); @assert length(methods(FS.prescribed_forcing_netcdf)) > 0; println("prescribed_forcing_netcdf symbol smoke test passed")'
    check_julia_smoke \
        "Julia smoke test: PlantHydraulicsModel API compatibility" \
        'import ClimaLand; import ClimaLand.Canopy as Canopy; import ClimaLand.Parameters as LP; import ClimaLand.Domains: Column; FT = Float64; domain = Column(; zlim = (FT(-9), FT(0)), nelements = 20, dz_tuple = (FT(1.5), FT(0.025)), longlat = (FT(11.64), FT(55.49))); canopy_domain = ClimaLand.Domains.obtain_surface_domain(domain); toml_dict = LP.create_toml_dict(FT); h = Canopy.PlantHydraulicsModel{FT}(canopy_domain, toml_dict); @assert hasproperty(h, :parameters); @assert !hasproperty(h, :n_stem); @assert !hasproperty(h, :n_leaf); println("PlantHydraulicsModel API smoke test passed")'
    check_julia_smoke \
        "Julia smoke test: AutotrophicRespiration ParamDict optional key handling" \
        'import ClimaLand; import ClimaLand.Canopy as Canopy; import ClimaLand.Parameters as LP; FT = Float64; toml_dict = LP.create_toml_dict(FT); p = Canopy.AutotrophicRespirationParameters(toml_dict); @assert p.Q10 > FT(0); println("AutotrophicRespiration ParamDict smoke test passed")'
    check_julia_smoke \
        "Julia smoke test: SoilCO2Model ParamDict optional key handling" \
        'import ClimaLand; import ClimaLand.Soil.Biogeochemistry as BGC; import ClimaLand.Parameters as LP; FT = Float64; toml_dict = LP.create_toml_dict(FT); p = BGC.SoilCO2ModelParameters(toml_dict); @assert p.Ea_sx > FT(0); @assert p.α_sx > FT(0); println("SoilCO2Model ParamDict smoke test passed")'
    check_julia_smoke \
        "Julia smoke test: SoilCO2 diffusivity and source compatibility" \
        'import ClimaLand; import ClimaLand.Soil.Biogeochemistry as BGC; import ClimaLand.Parameters as LP; FT = Float64; toml_dict = LP.create_toml_dict(FT); p = BGC.SoilCO2ModelParameters(toml_dict); Dco2 = BGC.co2_diffusivity(FT(298.15), FT(0.2), FT(101325), FT(0.35), FT(4.5), FT(0.45), p); Do2 = BGC.o2_diffusivity(FT(298.15), FT(0.2), FT(101325), FT(0.35), FT(4.5), FT(0.45), p); Rsm = BGC.microbe_source(FT(298.15), FT(0.2), FT(5.0), FT(0.1), p); @assert Dco2 > FT(0); @assert Do2 > FT(0); @assert Rsm >= FT(0); println("SoilCO2 diffusivity/source smoke test passed")'
    check_julia_smoke \
        "Julia smoke test: SoilCO2 Henry constant 3-arg compatibility" \
        'import ClimaLand; import ClimaLand.Soil.Biogeochemistry as BGC; FT = Float64; K = BGC.henry_constant(FT(1.0), FT(2000.0), FT(298.15)); @assert K > FT(0); println("SoilCO2 Henry constant smoke test passed")'
    check_julia_smoke \
        "Julia smoke test: SoilCO2 implicit aux compatibility" \
        'import ClimaLand; import ClimaLand.Soil.Biogeochemistry as BGC; @assert hasmethod(ClimaLand.make_update_implicit_aux, Tuple{BGC.SoilCO2Model}); println("SoilCO2 implicit aux smoke test passed")'
else
    fail "julia not found in PATH"
fi

echo
echo "[5/6] Checking SLURM dry-run parse (if available)"
if command -v sbatch >/dev/null 2>&1; then
    DRYRUN_OK=1
    for script in \
        experiments/calibrate_dk_sor/slurm_calibration.sh \
        experiments/callmip_uq_dk_sor/slurm_emulate_sample.sh \
        experiments/callmip_uq_dk_sor/slurm_callmip_simulations.sh \
        experiments/callmip_uq_dk_sor/slurm_postprocess.sh; do
        if sbatch --test-only "$script" >/tmp/callmip_sbatch_test.log 2>&1; then
            pass "sbatch --test-only: $script"
        else
            DRYRUN_OK=0
            fail "sbatch --test-only failed: $script"
            echo "See /tmp/callmip_sbatch_test.log"
        fi
    done
else
    fail "sbatch not found in PATH"
fi

echo
echo "[6/6] Summary"
echo "Pass: $PASS_COUNT"
echo "Fail: $FAIL_COUNT"

if [[ $FAIL_COUNT -eq 0 ]]; then
    echo "READY_TO_SUBMIT=YES"
    exit 0
else
    echo "READY_TO_SUBMIT=NO"
    exit 1
fi
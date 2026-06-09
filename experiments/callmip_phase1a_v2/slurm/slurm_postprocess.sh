#!/bin/bash
#SBATCH --job-name=dk_sor_postprocess
#SBATCH --output=experiments/callmip_phase1a_v2/slurm_logs/postprocess_%j.out
#SBATCH --error=experiments/callmip_phase1a_v2/slurm_logs/postprocess_%j.err
#SBATCH --partition=expansion
#SBATCH --account=esm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00

export CLIMACOMMS_DEVICE="CPU"
export CLIMACOMMS_CONTEXT="SINGLETON"

module use /groups/esm/modules
module load climacommon

REPO_ROOT="$SLURM_SUBMIT_DIR"
while [[ ! -d "$REPO_ROOT/.buildkite" && "$REPO_ROOT" != "/" ]]; do
    REPO_ROOT="$(dirname "$REPO_ROOT")"
done
cd "$REPO_ROOT"
echo "Running from: $PWD"

EXP=experiments/callmip_phase1a_v2
mkdir -p ${EXP}/slurm_logs
julia --project=${EXP} -e 'using Pkg; Pkg.instantiate()'
julia --project=${EXP} ${EXP}/analyze_posterior_ensemble.jl || echo "WARNING: posterior ensemble analysis failed"
julia --project=${EXP} ${EXP}/evaluate_calibration.jl       || echo "WARNING: evaluate_calibration failed"
julia --project=${EXP} ${EXP}/write_callmip_netcdf.jl
julia --project=${EXP} ${EXP}/plot_calibration_check.jl     || echo "WARNING: calibration plots failed"
julia --project=${EXP} ${EXP}/plot_posterior_distributions.jl || echo "WARNING: posterior plots failed"
echo "Postprocessing done at $(date)"

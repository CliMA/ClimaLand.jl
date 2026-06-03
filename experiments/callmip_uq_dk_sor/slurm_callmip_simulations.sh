#!/bin/bash
#SBATCH --job-name=callmip_sims
#SBATCH --partition=expansion
#SBATCH --account=esm
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/callmip_sims_%j.out
#SBATCH --error=logs/callmip_sims_%j.err

# =============================================================================
# CalLMIP Phase 1 – Prior + Posterior Simulations (DK-Sor)
#
# Runs two forward simulations (prior + posterior) via
# run_callmip_simulations.jl → callmip_model_interface.jl
#
# Outputs:
#   experiments/callmip_uq_dk_sor/output_callmip_sims/
#     iteration_000/member_001/callmip_diagnostics.jld2  (prior)
#     iteration_000/member_002/callmip_diagnostics.jld2  (posterior)
#
# After this job completes, run write_callmip_netcdf.jl to produce NetCDF.
# =============================================================================

set -euo pipefail

module use /groups/esm/modules
module load climacommon

export CLIMACOMMS_DEVICE="CPU"

# Resolve experiment directory from the sbatch submission directory
CLIMALAND_DIR="${SLURM_SUBMIT_DIR}"
SCRIPT_DIR="${SLURM_SUBMIT_DIR}/experiments/callmip_uq_dk_sor"

echo "ClimaLand directory : ${CLIMALAND_DIR}"
echo "Script directory    : ${SCRIPT_DIR}"
echo "SLURM job ID        : ${SLURM_JOB_ID}"
echo "SLURM nodes         : ${SLURM_NODELIST}"
echo "Start time          : $(date)"

# Create log directory
mkdir -p "${SCRIPT_DIR}/logs"

cd "${CLIMALAND_DIR}"

julia \
    --project=.buildkite \
    --threads="${SLURM_CPUS_PER_TASK:-8}" \
    -p 2 \
    experiments/callmip_uq_dk_sor/run_callmip_simulations.jl

echo "CalLMIP simulations finished at $(date)"
echo ""
echo "Next step — write CalLMIP NetCDF output:"
echo "  julia --project=experiments/callmip_uq_dk_sor \\"
echo "        experiments/callmip_uq_dk_sor/write_callmip_netcdf.jl"

#!/bin/bash
#SBATCH --job-name=dk_sor_gen_obs
#SBATCH --partition=expansion
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=dk_sor_gen_obs_%j.out
#SBATCH --error=dk_sor_gen_obs_%j.err

export CLIMACOMMS_DEVICE="CPU"
export CLIMACOMMS_CONTEXT="SINGLETON"

# Use ESM group modules (provides Julia 1.12.5 via climacommon)
module use /groups/esm/modules
module load climacommon

cd /home/renatob/ClimaLand.jl

julia --project=.buildkite -e 'using Pkg; Pkg.instantiate()'
julia --project=.buildkite experiments/calibrate_dk_sor/generate_observations.jl

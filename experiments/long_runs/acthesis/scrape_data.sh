#!/bin/bash
#SBATCH --job-name=acharbon_scrape_data
#SBATCH --output=data_scrape_output.txt
#SBATCH --error=data_scrape_error.txt
#SBATCH --time=72:00:00
#SBATCH --mail-user=acharbon@caltech.edu

python ClimaLand.jl/experiments/long_runs/acthesis/get_data.py
#!/bin/bash
#SBATCH --job-name=joco_job_test      # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=josephineconnelly@qgg.au.dk     # Where to send mail	
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=1gb                     # Job memory request
#SBATCH --time=00:05:00               # Time limit hrs:min:sec
#SBATCH --output=serial_test_%j.log   # Standard output 
#SBATCH --error serial_test_error_%j.log   # Standard error log
pwd; hostname; date

echo "Running R script on a single CPU core"

module load R

module load anaconda3/2020.11
$ conda create --name mshpr-env --channel conda-forge r-dplyr r-rmapshaper r-sf
$ conda activate mshpr-env
$ R

R devtools::install_github("DavisVaughan/furrr")



module load python

echo "Running plot script on a single CPU core"

python /data/training/SLURM/plot_template.py

date



wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
chmod +x miniconda.sh
bash miniconda.sh -b
./miniconda3/bin/conda init bash

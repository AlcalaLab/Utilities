#!/bin/bash
#SBATCH --job-name=sra_dta
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH -t 48:00:00
#SBATCH --mem=64G

python sra_download_trim_assemble.py sra_taxon_list.csv

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=500GB
#SBATCH --time=72:00:00
#SBATCH --partition=celltypes
source /allen/programs/celltypes/workgroups/hct/YiDing/miniconda3/bin/activate /home/giuseppe.saldi/miniconda3/envs/mvi
cd $SLURM_SUBMIT_DIR

ipython3 2_wrapatac.py

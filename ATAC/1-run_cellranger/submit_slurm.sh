#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH -p skylake
#SBATCH -A gottgens-sl2-cpu
#SBATCH --time=0-20:59:59
#SBATCH -e cellranger-err-%j.%N.err
#SBATCH --output cellranger-log-%J.txt

# -----------------------------------------------------------
#title           :cellranger_atac_count.sh
#author          :hpb29
#date            :20180404
#version         :0.7    

#description     :This script submits 10x libraries via
#                 SLURM to 'cellranger count' and expects 
#                 to be fed four positional arguments:
#                 $1 - library/sample id
#                 $2 - path to 10x reference
#                 $3 - path to fastq files dir
#                 $4 - expected number of cells

#usage            :On this file you are only supposed to change
#                 the SLURM parameters (if necessary). For the
#                 cellranger parameters it is best to submit
#                 through the 'run_cellranger.sh' file.
#                 
# ------------------------------------------------------------

echo "Started processing library $1 on `date`"
echo -e "JobID: $JOBID\n======"

cellranger-atac count --id=$1 \
                 --sample=$4 \
                 --reference=$2 \
                 --fastqs=$3

echo "Finishing processing library $1 on `date`"

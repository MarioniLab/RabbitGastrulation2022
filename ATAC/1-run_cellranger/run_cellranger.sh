#!/usr/local/bin/bash

# Example Usage
# ./run_cellranger -c lib/cellranger-atac-2.0.0/cellranger-atac-2.0.0 

while getopts : flag
do
    case "${flag}" in    
        r) ref=${OPTARG};;
		f) fastqdir=${OPTARG};;
		x) slurm=${OPTARG};;
    esac
done


# Run cellranger-atac on each sample
sbatch ${slurm} BGRGP1 ${ref} ${fastqdir} SLX18911_SINAA3,SLX19158_SINAA3 
sbatch ${slurm} BGRGP2 ${ref} ${fastqdir} SLX18911_SINAB3,SLX19158_SINAB3
sbatch ${slurm} BGRGP3 ${ref} ${fastqdir} SLX18911_SINAC3,SLX19158_SINAC3  
sbatch ${slurm} BGRGP4 ${ref} ${fastqdir} SLX18911_SINAD3,SLX19158_SINAD3 
sbatch ${slurm} BGRGP5 ${ref} ${fastqdir} SLX18911_SINAE3,SLX19158_SINAE3  
sbatch ${slurm} BGRGP6 ${ref} ${fastqdir} SLX18911_SINAF3,SLX19158_SINAF3
sbatch ${slurm} BGRGP7 ${ref} ${fastqdir} SLX18911_SINAG3,SLX19158_SINAG3
sbatch ${slurm} BGRGP8 ${ref} ${fastqdir} SLX18911_SINAH3,SLX19158_SINAH3
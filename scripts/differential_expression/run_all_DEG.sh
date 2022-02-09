#!/usr/local/bin/bash

# conda activate basic_renv
# cd /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/pipeline/evaluating
# ./run_all_DEG.sh -s /nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/pipeline/evaluating/run_DEG.R -f '/nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/gastruloids_scRNAseq/pipeline/processing/DEG_functions.R' -d '/nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/sce.RDS' -o '/nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/DE_out/' -g '/nfs/research/marioni/Leah/gastrulation_epigenetics/gastruloid_characterisation/data_dir/processed_files_Leah/scRNAseq/gene_metadata.txt.gz'

while getopts s:d:f:o:g: flag
do
    case "${flag}" in
        s) script=${OPTARG};;
        d) sce=${OPTARG};;
        f) functions=${OPTARG};;
        o) outbase=${OPTARG};;
        g) genemeta=${OPTARG};;
    esac
done

celltypes=('All' 'Mature Somites 2' 'NMP-derived Meso' 'Early Head Mesoderm' 'NMP' 'Spinal Cord' 'Pharyngeal Mesoderm?' 'Early NMP-independent Neural' 'Early NMP-independent Meso' 'Mature Endoderm' 'Late PS' 'Mid NMP-independent Meso' 'Late NMP-independent Neural' 'Late NMP-independent Meso' 'Caudal Epiblast' 'Inactive PS' 'NMP-derived Neural' 'Early PS' 'Early Mesp1 Mesoderm' 'Late Mesp1 Mesoderm' 'Endothelium' 'Separate FMH' 'Mature Somites 1' 'Late Head Mesoderm' 'Early Endoderm')
for ia in "${!celltypes[@]}"
do
    for ib in "${!celltypes[@]}"
    do
        if [ "$ia" -lt "$ib" ] ; then
            bsub -M 100000 "$script" -f "$functions" -s "$sce" -o "$outbase" -g "$genemeta" -a "${celltypes[$ia]}" -b "${celltypes[$ib]}"
            #echo "${celltypes[$ia]}" "${celltypes[$ib]}"
        fi
    done
done

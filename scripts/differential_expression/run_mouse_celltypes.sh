#!/usr/local/bin/bash

#./run_mouse_celltypes.sh -e r_v4 -s run_edgeR.R -d /nfs/research/marioni/dkeitley/RabbitGastrulation2021/data-out/compare_genes/m_downsample.rds -o /nfs/research/marioni/dkeitley/RabbitGastrulation2021/data-out/compare_genes/edgeR_results/mouse/ -l /nfs/research/marioni/dkeitley/RabbitGastrulation2021/data-out/compare_genes/edgeR_results/mouse/logs/

while getopts s:e:d:o:l: flag
do
    case "${flag}" in    
        s) script=${OPTARG};;
	e) env=${OPTARG};;
        d) sce=${OPTARG};;
        o) outdir=${OPTARG};;
	l) logdir=${OPTARG};;
    esac
done


celltypes=("Epiblast" "Primitive Streak" "Nascent mesoderm" "Haematoendothelial progenitors" "Anterior Primitive Streak" "PGC" "Allantois endothelium" "Visceral endoderm" "Cardiopharyngeal progenitors" "Embryo proper endothelium" "EMP" "Placodal ectoderm" "Lateral plate mesoderm" "Notochord" "NMPs" "Non-neural ectoderm 4" "Pharyngeal endoderm" "Dermomyotome" "MEP" "Cranial mesoderm" "Caudal mesoderm" "Mesenchyme" "Anterior somitic tissues" "Anterior cardiopharyngeal progenitors" "Cardiomyocytes FHF 1" "Venous endothelium" "Erythroid" "Mesothelium" "NMPs/Mesoderm-biased" "Cardiomyocytes FHF 2" "YS endothelium" "Sclerotome" "Presomitic mesoderm" "Allantois" "Endocardium" "Somitic mesoderm" "Epicardium" "Posterior somitic tissues" "Gut tube" "Nephron progenitors" "Non-neural ectoderm 2" "Limb mesoderm" "Thyroid primordium" "Migratory neural crest" "Megakariocytes" "Mesothelium-endothelium/Masked")

# Test of subset of cell types
#celltypes=("${celltypes[@]:0:5}")

for ia in ${!celltypes[@]}
do
    #echo ${celltypes[$ia]} ${celltypes[$ib]}

    celltypeA=${celltypes[$ia]////_}
    celltypeB="All"
    outstem=${celltypeA}_vs_${celltypeB}
    outstem=${outstem//[[:blank:]]/_}
    #echo $outstem

    stdout=${logdir}/${outstem}_stdout.txt
    stderr=${logdir}/${outstem}_stderr
    outfile=${outdir}/${outstem}.rds
    #echo $stdout

    bsub -q bigmem -M 128000 -o $stdout -e $stderr conda run -n $env Rscript "$script" -s "$sce" -o "$outfile" -g celltype -a "${celltypes[$ia]}" -b "All" -k sample

done

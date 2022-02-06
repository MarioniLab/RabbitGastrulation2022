#!/usr/local/bin/bash

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

mapfile -t celltypes < /nfs/research/marioni/dkeitley/RabbitGastrulation2021/data-in/compare_genes/rabbit_celltypes.csv

for ia in ${!celltypes[@]}
do
    for ib in ${!celltypes[@]}
    do
        if [ "$ia" -lt "$ib" ] ; then
	    stdout=$logdir/${celltypes[$ia]}_vs_${celltypes[$ib]}_stdout.txt
	    stderr=$logdir/${celltypes[$ia]}_vs_${celltypes[$ib]}_stderr
	    outfile=$outdir/${celltypes[$ia]}_vs_${celltypes[$ib]}.rds
            
	    bsub -M 100000 -o $stdout -e $stderr conda run -n $env Rscript "$script" -s "$sce" -o "$outfile" -g celltype -a "${celltypes[$ia]}" -b "${celltypes[$ib]}" -k sample
        fi
    done
done

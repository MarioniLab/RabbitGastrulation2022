#!/usr/local/bin/bash

#./run_nhood_tests.sh -e r_v4 -s nhood_tests.R -n /nfs/research/marioni/dkeitley/RabbitGastrulation2021/data-out/compare_nhoods/r_milo.rds -m /nfs/research/marioni/dkeitley/RabbitGastrulation2021/data-out/compare_nhoods/m_milo.rds -o /nfs/research/marioni/dkeitley/RabbitGastrulation2021/data-out/compare_nhoods/nhood_tests/ -l /nfs/research/marioni/dkeitley/RabbitGastrulation2021/data-out/compare_nhoods/nhood_tests/logs/

while getopts e:s:n:m:x:o:l: flag
do
    case "${flag}" in    
        e) env=${OPTARG};;
        s) script=${OPTARG};;
      	n) milo1=${OPTARG};;
	m) milo2=${OPTARG};;
	x) orthologs=${OPTARG};;
        o) outdir=${OPTARG};;
	l) logdir=${OPTARG};;
    esac
done


TESTS=("TEST_N_GENES" "TEST_GSPEC" "TEST_CORTYPE" "TEST_HVG_JOIN")

# Run a subset of tests
#TESTS=("${TESTS[@]:0:2}")

for ia in ${!TESTS[@]}
do

    test=${TESTS[$ia]}

    stdout=${test}_stdout.txt
    stderr=${test}_stderr

    bsub -q bigmem -M 256000 -o $stdout -e $stderr conda run -n $env Rscript "$script" -n "$milo1" -m "$milo2" -x $orthologs -t $test -o "$outdir" 

done

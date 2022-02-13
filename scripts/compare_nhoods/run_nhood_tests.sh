#!/usr/local/bin/bash

#./run_nhood_tests.sh -e r_v4 -s nhood_tests.R -n1 /nfs/research/marioni/dkeitley/RabbitGastrulation2021/data-out/compare_nhoods/r_milo.rds -n2 /nfs/research/marioni/dkeitley/RabbitGastrulation2021/data-out/compare_nhoods/m_milo.rds -o /nfs/research/marioni/dkeitley/RabbitGastrulation2021/data-out/compare_nhoods/tests/ -l /nfs/research/marioni/dkeitley/RabbitGastrulation2021/data-out/compare_nhoods/tests/logs/

while getopts e:s:n1:n2:o:l: flag
do
    case "${flag}" in    
        e) env=${OPTARG};;
        s) script=${OPTARG};;
      	n1) milo1=${OPTARG};;
		n2) milo2=${OPTARG};;
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

    bsub -q bigmem -M 256000 -o $stdout -e $stderr conda run -n $env Rscript "$script" -n1 "$milo1" -n2 "$milo2" -t $test -o "$outdir" 

done

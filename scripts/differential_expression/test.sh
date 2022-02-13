#!/usr/local/bin/bash

mapfile -t my_array < /nfs/research/marioni/dkeitley/RabbitGastrulation2021/data-in/compare_genes/rabbit_celltypes.csv
printf '%s\n' "${my_array[@]}"

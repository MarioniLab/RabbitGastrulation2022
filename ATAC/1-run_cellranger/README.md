# Run cellranger-ATAC

This directory contains the scripts used to run the 10X Genomics cellranger-ATAC pipeline on the scATAC-seq libraries produced for the rabbit.


## Requirements
The cellranger-ATAC software can be downloaded from the 10X Genomics website [here](https://support.10xgenomics.com/single-cell-atac/software/downloads/latest) and installed using instructions [here](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/installation). 


The sequencing files can be downloaded from ArrayExpress using accession [E-MTAB-11804](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11804/). 

The cellranger-ATAC pipeline requires constructing a reference package from a reference genome sequence and gene annotations. The rabbit reference genome an gene annotations used can be accessed [here](https://content.cruk.cam.ac.uk/jmlab/RabbitGastrulation2022/data/). 

The 10X reference package can then be created using 
```

```



 ## Usage
 
The `run_cellranger.sh` script submits jobs via [SLURM](https://slurm.schedmd.com/overview.html) to execute the `cellranger-atac count` pipeline on scATAC-seq samples.

Modifications are required in order to be compatible with other cluster scheduling systems. 

The SLURM script `submit_slurm.sh` can be edited to specify different cluster settings.

Input files are provided to `run_cellranger.sh` using the following flags:

```
-r	The directory of the reference
-f 	The directory containing ATAC-seq fastq files
-x 	The path to the SLURM submission script (e.g. submit_slurm.sh)
```

The `submit_slurm.sh` script assumes that cellranger-ATAC is available on the PATH environment variable. E.g. 

```
export PATH=$SCRABBIT_HOME/lib/cellranger-atac-2.0.0/cellranger-atac-2.0.0:$PATH 
```


### Example
```
./run_cellranger.sh -r $SCRABBIT_HOME/ATAC/data-in/run_cellranger/10X-rabbit-reference/ -f $SCRABBIT_HOME/ATAC/data-in/run_cellranger/fastqs/ -x submit_slurm.sh
```
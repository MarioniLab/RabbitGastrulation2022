# Create ArchR project

This directory contains scripts to create an [ArchR](https://www.archrproject.com/) project from the raw cellranger-ATAC output files. 


## Requirements
We ran these scripts using [Snakemake](https://snakemake.readthedocs.io/en/stable/), which should be installed as part of the `scrabbit-atac` conda environment (see the project's `envs/` directory). Alternatively, installation instructions can be found [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

The scripts can also be run independently outside of Snakemake. Feel free to file an issue if you have problems with this. 

The pipeline was ran in a cluster environment using [SLURM](https://slurm.schedmd.com/).   The `run_snakemake.sh` file will need to be modified in order to work with other scheduling systems.

The scripts rely on R packages that are installed with [scrabbitr](https://github.com/dkeitley/scrabbitr). 


## Usage
```
./run_snakemake.sh
```
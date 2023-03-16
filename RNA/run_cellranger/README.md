# Run Cell Ranger

This directory contains scripts used to run the 10X Genomics Cell Ranger pipleine on the scRNA-seq libraries produced from the rabbit samples.

## Requirements

Cell Ranger can be downloaded and installed from the 10X Genomics website - [link](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest).

The reference genome can be downloaded from [Ensembl](https://www.ensembl.org/Oryctolagus_cuniculus/Info/Index).

Our modified GTF file can be downloaded from [here](https://content.cruk.cam.ac.uk/jmlab/RabbitGastrulation2022/data/RNA/genome_annotation/).

The scRNA-seq raw fastq files can be downloaded from the ArrayExpress website using accession number [**E-MTAB-11836**](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11836/).

## Usage

The first step is to create a Cell Ranger reference package using the reference genome and gene annotations. The code for doing this is available in `create_reference.md`

The `cellranger count` pipeline is then run for each sample providing an id, the reference package created above, the raw fastq files, sample IDs and the number of cells expected from the experiment.

e.g.

    $SCRABBIT_HOME/lib/cellranger-3.1.0/cellranger count --id=SLX18995_rabbit_SIGAA9 
    --transcriptome=data-out/rabbit-cr-reference/filtered/ 
    --fastqs=$SCRABBOT_HOME/data-in/rabbit/fastqs 
    --sample=SIGAA9,B2_SIGAA9,SLX18995_B3/B3_SIGAA9,SLX18995_B4/SIGAA9 
    --expect-cells=4174

See Cell Ranger [documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count#cr-count) for more info.

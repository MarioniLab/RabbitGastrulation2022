# Improving gene annotations

These scripts were designed to investigate the alignment of 10X sequencing reads to the rabbit reference genome and to the modify the genome annotation.

[annotate.R](src/annotate.R) - Used to visualise distances of reads to their nearest annotated gene. Writes artificial GTF annotation entries extending exons from the 3' end.

[get_alignments.pl](src/get_alignments.pl) - Queries Ensembl Compara for human/rabbit alignments. Exports alignments of 3' exons and UTRs which can be used to suggest regions to annotate.

## Requirements

The OryCun2.0 reference genome can be obtained from [Ensembl](https://www.ensembl.org/Oryctolagus_cuniculus/Info/Index).

`annotate.R` visualises intergenic reads obtained from the cellranger BAM output file for the SIGAC11 GD8 sample. Intergenic reads were obtained by filtering on the `'RE'` BAM alignment tag using [SAMtools](http://www.htslib.org/) (see details [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam#bam-align-tags)).

    (samtools view -H possorted_genome_bam.bam; samtools view possorted_genome_bam.bam | grep -w 'RE:A:I') | samtools view -bS - > intergenic_reads.bam

This `get_alignments.pl` script utilises the Ensembl Perl API. Please refer to the following links for installation instructions and tutorials.

[Installation instructions for the Perl API](http://www.ensembl.org/info/docs/api/api_installation.html)

[Getting started with the Perl API](http://www.ensembl.org/info/docs/api/general_instructions.html)

[Accessing Ensembl Compara with the Perl API](http://www.ensembl.org/info/docs/api/compara/compara_tutorial.html)

[Ensembl Compara API Documentation](http://www.ensembl.org/info/docs/Doxygen/compara-api/index.html)

Alignments were obtained for genes which have a one-to-one orthology relationship between the mouse and human but which were missing from the rabbit reference. This list is provided in `data-in/hs_mm_candidates.tsv` and was obtained using Ensembl [biomaRt](<https://bioconductor.org/packages/release/bioc/html/biomaRt.html>).

E.g.

    # Get mouse-human homolog data
    mm_hs_homologs <- scrabbitr::getEnsemblHomologs("mmusculus", mouse_genes$ensembl_id,"hsapiens")

    # Filter mouse-human homologs
    mm_hs_homologs <- mm_hs_homologs[mm_hs_homologs[["hsapiens_homolog_ensembl_gene"]]!="",]

    # Get mouse-rabbit homolog data
    mm_hs_genes <- unique(mm_hs_homologs$ensembl_gene_id)
    mm_oc_homologs <- getEnsemblHomologs("mmusculus", mm_hs_genes,"ocuniculus")

    # Filter mouse-human not rabbit homologs
    mm_hs_not_oc <- mm_oc_homologs[mm_oc_homologs[["ocuniculus_homolog_ensembl_gene"]]=="",]

    # Export
    out <- mm_hs_homologs[mm_hs_homologs$ensembl_gene_id %in% mm_hs_not_oc$ensembl_gene_id,]
    write.table(out, file="data-out/hs_mm_candidates.tsv", quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)

### Background

In the reverse transcription step of the 3' 10X single-cell RNA-sequencing protocol, primers containing cellular and molecular barcodes hybridise to the poly-A tails at the 3' end of mRNA molecules. cDNA synthesis then occurs from this 3' poly-A tail in the 3'-5' direction, capturing the adjacent \~98 base pairs of the transcript sequence.

Often in non-model organisms, these 3' ends of transcripts are less confidently annotated. This means that large numbers of 10X sequencing reads may align to unannotated regions of the reference genome, just a short distance away from the 3' ends of annotated genes. These reads will not be counted by the 10X cellranger pipeline.

[FOXA2 annotation](res/foxa2_annotation.jpg)

The workaround used here, is to artificially extend the 3' annotations of genes by a short distance so that these reads now overlap annotated regions and thus are counted by the cellranger pipeline. Functions to do this are provided in the [annotate.R](src/annotate.R) script.

A second challenge is that many well-known genes are missing from the OryCun2.0 reference. One example is Six3, an anterior marker that is well-conserved across bilaterian species. Here, we have used alignment data from Ensembl Compara to add annotations of genes that may bey present in the genome but which are missing from the reference annotation. The [get_alignments.pl](src/get_alignments.pl) script provides methods to query Ensembl rabbit-human alignment data to provide candidates for new annotations.

[FOXA2 annotation](res/six3_annotation.jpg)

See the paper methods for more details.

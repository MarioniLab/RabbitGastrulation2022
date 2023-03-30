# Make genome reference

Filtering GTF files prior to compiling reference genome for O. cuniculus.

    awk '$1 !~ /^#/' extension_plus_alignments.gtf | cut -f 1 | sort | uniq | screen_list.pl - Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa keep > filtered_01042020_Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa

Now we count the lines in the two fasta files, we see that we filtered out over half of the contigs, however the file size is about the same in both cases, 2.6 GB:

    grep ">" filtered_01042020_Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa | wc -l
    #1426

    grep ">" Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa | wc -l
    #3242

Next we filter the GTF (see [10X website](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#mkgtf)).

    /rds/project/bg200/rds-bg200-hphi-gottgens/users/mlnt2/cellranger/cellranger-3.1.0/cellranger mkgtf \
    extension_plus_alignments.gtf \
    Oryctolagus_cuniculus.OryCun2.0.99.extension_plus_alignments_filtered.gtf \
                       --attribute=gene_biotype:protein_coding \
                       --attribute=gene_biotype:lincRNA \
                       --attribute=gene_biotype:antisense \
                       --attribute=gene_biotype:IG_LV_gene \
                       --attribute=gene_biotype:IG_V_gene \
                       --attribute=gene_biotype:IG_V_pseudogene \
                       --attribute=gene_biotype:IG_D_gene \
                       --attribute=gene_biotype:IG_J_gene \
                       --attribute=gene_biotype:IG_J_pseudogene \
                       --attribute=gene_biotype:IG_C_gene \
                       --attribute=gene_biotype:IG_C_pseudogene \
                       --attribute=gene_biotype:TR_V_gene \
                       --attribute=gene_biotype:TR_V_pseudogene \
                       --attribute=gene_biotype:TR_D_gene \
                       --attribute=gene_biotype:TR_J_gene \
                       --attribute=gene_biotype:TR_J_pseudogene \
                       --attribute=gene_biotype:TR_C_gene

    /rds/project/bg200/rds-bg200-hphi-gottgens/users/mlnt2/cellranger/cellranger-3.1.0/cellranger mkref \
    --genome=rabbit_filtered \
    --fasta=filtered_01042020_Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa \
    --genes=Oryctolagus_cuniculus.OryCun2.0.99.extension_plus_alignments_filtered.gtf

Reference genome stored at: `"/rds/project/bg200/rds-bg200-hphi-gottgens/users/mlnt2/rabbit_reference/rabbit_filtered/"`

## ATAC

Check out the [10X ATAC website](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/advanced/references#mkref) on how to create custom reference file for rabbit ATACseq.

    /rds/project/bg200/rds-bg200-hphi-gottgens/users/mlnt2/cellranger/cellranger-atac-1.2.0/cellranger-atac mkref ocun_atac --config=ocun.config

    {
        GENOME_FASTA_INPUT: "/rds/project/bg200/rds-bg200-hphi-gottgens/users/mlnt2/rabbit_reference/filtered_Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa",
        GENE_ANNOTATION_INPUT: "/rds/project/bg200/rds-bg200-hphi-gottgens/users/mlnt2/rabbit_reference/Oryctolagus_cuniculus.OryCun2.0.99.filtered.gtf",
        MOTIF_INPUT: "",
        ORGANISM: "Oryctolagus cuniculus",
        PRIMARY_CONTIGS: ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "X"],
        NON_NUCLEAR_CONTIGS: ["MT"]
    }

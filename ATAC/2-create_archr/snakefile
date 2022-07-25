import os
from re import search
import getpass


############
## Config ##
############
host = os.uname()[1]
configfile: "/rds/project/rds-SDzz0CATGms/users/bt392/04_rabbit_scATAC/code/config.yaml"


###########
## Rules ##
###########

rule all:
    input:
        expand(config["directories"]["processed_data"]+"/rabbit_{sample}.arrow", sample=config["samples"]), # Create Arrows
        config["directories"]["processed_data"]+"/02_completed.txt", # Create ArchR project  
        #config["directories"]["results"]+"/qc/sample_metadata_after_qc.txt.gz", # QC  
        #expand("%s/dimensionality_reduction/lsi_{matrix}_nfeatures{nfeatures}_ndims{ndims}.txt.gz" % (config["directories"]["results"]),
        #    dimred_matrix = config["dimensionality_reduction"]["dimred_matrix"], 
        #    dimred_nfeatures = config["dimensionality_reduction"]["dimred_nfeatures"], 
        #    dimred_ndims = config["dimensionality_reduction"]["dimred_ndims"]
        #    ) # Tile matrix dim red first pass
        #expand(config["directories"]["results"]+"/dimensionality_reduction/{dimred_matrix}_final.txt", matrix= config["dimensionality_reduction_Tile"]["matrix"])
        

########################
## Create arrow files ##
########################

rule create_arrow_files:
    input:
        script = config["scripts"]["create_arrow_files"],
        fragments_files = config["directories"]["original_data"]+"/{sample}_fragments.tsv.gz"
    output:
        config["directories"]["processed_data"]+"/rabbit_{sample}.arrow"
    params:
#        outdir = config["directories"]["processed_data"],
        sample = config["samples"],
        min_fragments = config["create_arrow_files"]["min_fragments"],
        min_tss_score = config["create_arrow_files"]["min_tss_score"]
    threads: 
        config["slurm"]["create_arrow_files"]["threads"]
    resources:
        mem_mb = config["slurm"]["create_arrow_files"]["memory"]
    log: 
        "logs/create_arrow_files_{sample}.log"
    shell:
        "Rscript {input.script}  --sample {wildcards.sample}  --min_fragments {params.min_fragments} --min_tss_score {params.min_tss_score}  > {log}" # --fragments_files {input.fragments_files} --outdir {params.outdir} --threads {threads} 


##########################
## Create ArchR project ##
##########################

rule create_archr_project:
    input:
        script = config["scripts"]["create_archr_project"],
        arrow_files = expand(rules.create_arrow_files.output, sample=config["samples"])
    output:
        config["directories"]["processed_data"]+"/02_completed.txt"
    params:
        outdir = config["directories"]["processed_data"]
    threads: 
        config["slurm"]["create_archr_project"]["threads"]
    resources:
        mem_mb = config["slurm"]["create_archr_project"]["memory"]
    log: 
        "logs/create_archr_project.log"
    shell:
        "Rscript {input.script} --arrow_files {input.arrow_files} > {log}"

########
## QC ##
########

rule qc_archr:
    input:
        script = config["scripts"]["qc_archr"],
        proj = rules.create_archr_project.output
    output:
        config["directories"]["results"]+"/qc/qc_FragmentSizeDistribution.txt.gz",
        config["directories"]["results"]+"/qc/qc_FragmentSizeDistribution.pdf",
        config["directories"]["results"]+"/qc/qc_TSSenrichment.txt.gz",
        config["directories"]["results"]+"/qc/qc_TSSenrichment.pdf",
        config["directories"]["results"]+"/qc/qc_metrics_histogram.pdf",
        config["directories"]["results"]+"/qc/qc_metrics_barplot.pdf",
        metadata=config["directories"]["results"]+"/qc/sample_metadata_after_qc.txt.gz"
    params:
        min_tss_score = config["qc_archr"]["min_tss_score"],
        min_log_nFrags = config["qc_archr"]["min_log_nFrags"]
    threads: 
        config["slurm"]["qc_archr"]["threads"]
    resources:
        mem_mb = config["slurm"]["qc_archr"]["memory"]
    log: 
        "logs/qc_archr.log"
    shell:
        "Rscript {input.script} --min_tss_score {params.min_tss_score} --min_log_nFrags {params.min_log_nFrags} > {log}"
        
        
##############################
## Dimensionality reduction ##
##############################

#rule dimensionality_reduction_Tile: 
#    input:
#        script = config["scripts"]["dimensionality_reduction"]
#    output:
#        token = config["directories"]["results"]+"/dimensionality_reduction/TileMatrix_first_pass.txt"
#    params:
#        seed = 42,
#        outdir = config["directories"]["results"] + "/dimensionality_reduction",
#        matrix = config["dimensionality_reduction_Tile"]["matrix"],
#        nfeatures = config["dimensionality_reduction_Tile"]["nfeatures"],
#        ndims = config["dimensionality_reduction_Tile"]["ndims"],
#        n_neighbors = config["dimensionality_reduction_Tile"]["n_neighbors"],
#        min_dist = config["dimensionality_reduction_Tile"]["min_dist"],
#        colour_by = config["dimensionality_reduction_Tile"]["colour_by"]
#    threads: 
#        config["slurm"]["dimensionality_reduction_Tile"]["threads"]
#    resources:
#        mem_mb = config["slurm"]["dimensionality_reduction_Tile"]["memory"]
#    log: 
#        "logs/umap_{dimred_matrix}_nfeatures{dimred_nfeatures}_ndims{dimred_ndims}.log"
#    shell:
#        "Rscript {input.script} --matrix {params.matrix} --nfeatures {wildcards.nfeatures} --ndims {wildcards.ndims} \
#        --n_neighbors {wildcards.n_neighbors} --min_dist {wildcards.min_dist} --colour_by {params.colour_by} \
#        --seed {params.seed} --outdir {params.outdir} > {log}"

#rule dimensionality_reduction_Tile_final: 
#    input:
#        script = config["scripts"]["dimensionality_reduction"],
#        test_previous = rule.dimensionality_reduction_Tile.output
#    output:
#        token = config["directories"]["results"]+"/dimensionality_reduction/TileMatrix_final.txt"
#    params:
#        seed = 42,
#        outdir = config["directories"]["results"] + "/dimensionality_reduction",
#        matrix = "TileMatrix",
#        nfeatures = config["dimensionality_reduction_Tile_final"]["nfeatures"],
#        ndims = config["dimensionality_reduction_Tile_final"]["ndims"],
#        n_neighbors = config["dimensionality_reduction_Tile_final"]["n_neighbors"],
#        min_dist = config["dimensionality_reduction_Tile_final"]["min_dist"],
#        colour_by = config["dimensionality_reduction_Tile_final"]["colour_by"]
#   threads: 
#        config["slurm"]["dimensionality_reduction_Tile_final"]["threads"]
#    resources:
#        mem_mb = config["slurm"]["dimensionality_reduction_Tile_final"]["memory"]
#    log: 
#        "logs/umap_{dimred_matrix}_nfeatures{dimred_nfeatures}_ndims{dimred_ndims}.log"
#    shell:
#        "Rscript {input.script} --matrix {params.matrix} --nfeatures {wildcards.nfeatures} --ndims {wildcards.ndims} \
#        --n_neighbors {params.n_neighbors} --min_dist {params.min_dist} --colour_by {params.colour_by} \
#        --seed {params.seed} --outdir {params.outdir} > {log}"


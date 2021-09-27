
#' Checks if provided genes are in Ensembl gene ID format.
#'  @export
isEnsemblID <- function(genes) {
  return(grepl("ENS[A-Z]+[0123456789]{11}",x=genes))
}



#' Queries Ensembl BioMart database for homologs between species.
#'@export
getEnsemblHomologs <- function(ref_species,ref_genes,query_species) {

  if(!any(isEnsemblID(ref_genes))) {
    stop("Genes given are not valid Ensembl IDs")
  }

  print("Retrieving homologs from Ensembl. This may take a few minutes...")

  mart = useMart("ensembl", dataset = paste0(ref_species,"_gene_ensembl"))

  homolog_data = getBM(attributes = c("ensembl_gene_id", "external_gene_name","chromosome_name",
                                      "start_position","end_position","strand",
                                      paste0(query_species,"_homolog_ensembl_gene"),
                                      paste0(query_species,"_homolog_associated_gene_name"),
                                      paste0(query_species,"_homolog_chromosome"),
                                      paste0(query_species,"_homolog_chrom_start"),
                                      paste0(query_species,"_homolog_chrom_end"),
                                      paste0(query_species,"_homolog_orthology_type"),
                                      paste0(query_species,"_homolog_orthology_confidence"),
                                      paste0(query_species,"_homolog_perc_id"),
                                      paste0(query_species,"_homolog_perc_id_r1"),
                                      paste0(query_species,"_homolog_goc_score"),
                                      paste0(query_species,"_homolog_wga_coverage")),
                       filters = "ensembl_gene_id",
                       values = ref_genes,
                       mart = mart)

  return(homolog_data)

}



#' Obtain one-to-one orthologs from Ensembl.
#' @param query_genes Genes from which to find one-to-one orthologs.
#' @param query_species Species to find one-to-one orthologs with provided in
#' Ensembl format (e.g. "hsapiens", "mmusculus").
#' @param one2one.data A data.frame of results retrieved from the Ensembl BioMart
#' database containing "_homolog_ensembl_gene" or  "_homolog_associated_gene_name"
#' fields.
#' @param outputIDs Specifies whether to return the one-to-one orthologs in
#' Ensembl gene ID format or using their gene names.
#' @export
# one2one.data assumed to be in format of getHomologData
# query_species in Ensembl format (e.g. hsapiens, mmusculus)
filterSharedFeatures <- function(query_genes, query_species,one2one.data,outputIDs=TRUE) {
  id <- TRUE
  if(!any(isEnsemblID(query_genes))) {
    # Assume query_genes given as gene name
    id <- FALSE
    field <- paste0(query_species,"_homolog_associated_gene_name")
  } else {
    field <- paste0(query_species,"_homolog_ensembl_gene")
  }

  one2ones_in_dataset <- one2one.data[[field]] %in% query_genes
  query.one2ones <- one2one.data[[field]][one2ones_in_dataset]

  ref.one2ones <- one2one.data$ensembl_gene_id[one2ones_in_dataset]

  if(!outputIDs) {
    gene_names <- one2one.data$external_gene_name[one2ones_in_dataset]
    ref.one2ones[gene_names!=""] <- gene_names[gene_names!=""]
  }

  sfs <- data.frame(ref=sapply(ref.one2ones,as.character),
                    query=sapply(query.one2ones,as.character),stringsAsFactors = FALSE)

  return(sfs)

}

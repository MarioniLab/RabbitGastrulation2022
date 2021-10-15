

# TODO: Add to scrabbitr
getGeneMeans <- function(gene, mat, celltypes){
  expr <- mat[gene,]
  medians = aggregate(formula = expr ~ celltypes, FUN = mean)
  out = matrix(medians$expr, nrow = 1, dimnames = list(gene, medians$celltypes))
  return(out)
}


#' @export
getMeanExpression <- function(sce,genes,group_by="celltype",norm="min-max") {
  genes <- genes[!is.na(genes)]
  mat <- logcounts(sce)[genes,]
  means <- lapply(genes, getGeneMeans, mat = mat, celltypes = sce[[group_by]])
  combined <- do.call(rbind,means)
  celltypes <- colnames(combined)
  if(norm=="min-max") {
    combined <- sweep(combined, 1, apply(combined, 1, max), "/")
  } else if(norm=="z-scale") {
    combined <- t(apply(combined, 1, scale))
    colnames(combined) <- celltypes
  } 
  return(combined)
}




#' Get dataframe of shared and unique markers.
#' Used for mouse/rabbit cell type markers heatmap. 
#' @export
getSharedUniqueMarkers <- function(merge.markers,N=1) {
  
  rchosen <- NULL
  mchosen <- NULL
  
  ncelltypes <- length(names(merge.markers))
  rshared <- mshared <- runiq <- muniq <- setNames(data.frame(matrix(ncol = ncelltypes,
                                                                     nrow = 1)),
                                                   names(merge.markers))
  for (celltype in names(merge.markers)) {
    markers <- merge.markers[[celltype]]
    
    # Make sure genes not already chosen
    markers <- markers[!(markers$rgene%in%rchosen) &
                         !(markers$mgene%in%mchosen),]
    
    markers <- markers[order(pmax(markers[,"FDR.x"],
                                  markers[,"FDR.y"]),
                             decreasing=F),]
    
    rshared_markers <- markers[markers$class=="shared","rgene"][1:N]
    mshared_markers <- as.character(markers[markers$class=="shared","mgene"][1:N])
    
    markers <- markers[order(markers[,"FDR.x"],
                             decreasing=F),]
    runiq_markers <- markers[markers$class=="rabbit unique","rgene"][1:N]
    
    markers <- markers[order(markers[,"FDR.y"],
                             decreasing=F),]
    muniq_markers <- markers[markers$class=="mouse unique","mgene"][1:N]
    
    rshared[[celltype]] <- rshared_markers
    mshared[[celltype]] <- mshared_markers
    runiq[[celltype]] <- runiq_markers
    muniq[[celltype]] <- muniq_markers
    
    rchosen <- c(rchosen,rshared_markers,runiq_markers)
    mchosen <- c(mchosen,mshared_markers,muniq_markers)
  }
  
  return(list(rshared=rshared,mshared=mshared,runiq=runiq,muniq=muniq))
  
}



#' Used to obtain genes related to apoptosis and proliferation for Waddington-OT analysis
#' Positively regulates...
#'   cell_death <- "GO:0010942"
#'   cell_population_proliferation <- "GO:0008284"
#'   
#' @importFrom biomaRt getBM useMart
#' @export
getGenesFromGOTerms <- function(species,go_terms,return_value="ensembl_gene_id") {
  mart <- biomaRt::useMart("ensembl", dataset = paste0(species,"_gene_ensembl"))
  
  data = biomaRt::getBM(attributes=c("ensembl_gene_id", "external_gene_name","go_id","name_1006",
                            "namespace_1003",
                            "description","gene_biotype"),
               filters="go",
               values = go_terms,
               uniqueRows = TRUE,
               mart=mart)
  
  genes <- unique(data[[return_value]])
  
  return(genes)
  
}

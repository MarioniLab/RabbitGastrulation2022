

#' Uniformly downsample a SingleCellExperiment object
#'
#' @param sce SingleCellExperiment object to downsample
#' @param ncells Number of cells required
#' @return Downsampled SingleCellExperiment object
#' @export
downsampleSCE <- function(sce, ncells) {
  cells <- sample(seq(1,ncol(sce),by=1),ncells,replace=FALSE)
  return(sce[,cells])
}



#' Write colData to a GMT file
#'
#' Used to export cell sets for Waddington-OT analysis.
#'
#' @details The gmt format consists of one cell set per line. Each line is a
#' tab-separated list composed as follows : \cr
#' \itemize{
#' \item The set name (can contain spaces) \cr
#' \item A commentary / description of the set (may be empty or contain spaces) \cr
#' \item A tab-separated list of set members
#' }
#'
#'
#'
#' @param sce SingleCellExperiment object
#' @param group colData observation to export
#' @param gmt_path Path to write GMT file
#' @export
exportToGMT <- function(sce,group,gmt_path) {
  lapply(unique(sce[[group]]),function(x) {
    out <- data.frame(c(x, "-", colnames(sce)[sce[[group]] == x]))
    write.table(t(out), gmt_path, sep = "\t", col.names = F,
                append = T,row.names = F,quote=F)
  })
}







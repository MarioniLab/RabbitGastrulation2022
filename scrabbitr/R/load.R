

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



#' Load csv files exported from WOT trajectory/fate AnnData
#' @export
loadWOTData <- function(wot_path) {
  traj <- read.csv(paste0(wot_path,"X.csv"),
                   sep="\t",header=FALSE)

  traj_obs <- read.csv(paste0(wot_path,"obs.csv"),
                       sep="\t",header=TRUE)
  traj_var <- read.csv(paste0(wot_path,"var.csv"),
                       sep="\t",header=FALSE)

  rownames(traj) <- traj_obs$X
  colnames(traj) <- traj_var[,1]
  traj <- cbind(traj,day=traj_obs$day)

  return(traj)

}



#' Load macaque data from Yang et al. 2021
#' @export
loadYang2021 <- function(data_path) {

  # Load meta
  meta_raw <- readRDS(paste0(data_path,"metadata_full.Rds"))
  meta_df <- as.data.frame(dplyr::bind_rows(meta_raw))
  rownames(meta_df) <- meta_df$cell

  # Load counts
  counts <- fread(file=paste0(data_path,"CM_filtered_SoupX_corrected.tsv"),sep="\t",
                  data.table=F)
  genes<-counts$V1
  counts <- counts[,grepl("WT*",colnames(counts))]
  counts <- sapply(counts,FUN=as.numeric)
  rownames(counts) <- genes
  logcounts <- Matrix::Matrix(counts,sparse=TRUE)


  umap_df <- meta_df[colnames(counts),c("UMAP_1","UMAP_2")]
  meta <- meta_df[colnames(counts),c("cell","lineage")]
  colnames(meta)[2] <- "celltype"

  # Add stage info
  meta$stage <- ""
  meta$stage[grepl("WT_d10_*", meta$cell)] <- "d10"
  meta$stage[grepl("WT_d12_*", meta$cell)] <- "d12"
  meta$stage[grepl("WT_d14_*", meta$cell)] <- "d14"

  meta <- meta[,-1]

  sce <- SingleCellExperiment(assays=list(counts=counts),colData=meta,
                              reducedDims=list(UMAP=umap_df))
  return(sce)
}



#' Load human data from Tyser et al. 2021
#' @export
loadTyser2020 <- function(data_path) {
  X <- readRDS(paste0(data_path,"expression_values.rds"))
  X_mat <- t(as.matrix(X))
  colnames(X_mat) <- paste0("cell_",seq(1,ncol(X_mat)))
  meta <- readRDS(paste0(data_path,"annot_umap.rds"))
  umap <- as.data.frame(meta[,c("X0","X1")])

  sce <- SingleCellExperiment(assays=list(logcounts = Matrix(X_mat,sparse = TRUE)),
                              colData = meta,
                              reducedDims = list(UMAP = umap))
  sce$stage <- "16dpf"

  colnames(colData(sce))[colnames(meta) == "cluster_id"] <- "celltype"
  return(sce)
}




#' @export
loadRabbitData <- function(data_path, normalise=TRUE) {

  # Load counts
  sce <- readRDS(paste0(data_path,"normalised_full.RDS"))

  # Load row data
  genes <- read.csv(paste0(data_path,"features.tsv"),sep="\t",header = F,
                    col.names = c("ensembl_id","gene_name","feature_type"))

  # Load metadata
  meta <- read.csv(paste0(data_path,"meta.tab"),
                   sep = "\t",header = TRUE)
  rownames(meta) <- meta$index
  celltypes <- read.csv(paste0(data_path,"annotation_12-07-21.tsv"),sep="\t",
                        row.names = 1)
  sce$celltype <-  celltypes[colnames(sce),"updated_celltype"]

  sce$cluster <- meta[,"leiden_2"]


  # Load corrected PCs
  pcs <- data.table::fread(paste0(data_path,"corrected_pcs.tsv"),sep="\t",
               data.table = F,col.names=c("cell",paste0("PC",seq(1:50))),showProgress = F,fill=T)
  rownames(pcs) <- pcs$cell
  pcs <- pcs[,-1]
  reducedDim(sce,"PCA") <- pcs[colnames(sce),]


  # Load umap
  umap <- read.table(paste0(data_path,"umap.tsv"),sep="\t")
  rownames(umap) <- meta$index
  reducedDim(sce,"UMAP") <- umap[colnames(sce),]

  # Apply normalisation
  if(normalise) {
    size_factors <- read.csv(paste0(data_path,"size_factors.tsv"),
                             sep="\t",header = T,  col.names= c("cell","sizeFactors"),row.names = "cell" )
    sce <- scater::logNormCounts(sce,size_factors=size_factors[colnames(sce),"sizeFactors"])

  }

  return(sce)

}

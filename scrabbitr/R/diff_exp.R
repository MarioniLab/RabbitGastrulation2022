

#' Prepare groupings in SingleCellExperiment object for running edgeR
#' @export
prepareEdgeR <- function(sce, group_by, group1, group2, block_by) {

  # Prepare groups
  if ((group1 != 'All') & (group2 != 'All')) {
    sce_sub <- sce[,colData(sce)[[group_by]] %in% c(group1, group2)]
    sce_sub$group <- colData(sce_sub)[[group_by]]
  } else if (group1 == 'All') {
    sce_sub <- sce
    sce_sub$group <- as.character(colData(sce_sub)[[group_by]])
    sce_sub$group[sce_sub$group != group2] <- 'All'
  } else if (group2 == 'All') {
    sce_sub <- sce
    sce_sub$group <- as.character(colData(sce_sub)[[group_by]])
    sce_sub$group[sce_sub$group != group1] <- 'All'
  }

  sce_sub$batch <- colData(sce_sub)[[block_by]]
  groups <- unique(sce_sub$group)
  sce_sub$group <- factor(sce_sub$group, levels = groups)

  return(sce_sub)
}


#' Compute differentially expressing genes using EdgeR
#' @importFrom data.table data.table setnames
#' @importFrom edgeR glmQLFit estimateDisp topTags glmQLFTest
#' @importFrom dplyr %>%
#' @importFrom DelayedArray rowMeans
#' @export
runEdgeR <- function(sce, min_detection_rate_per_group=0.1, min.logFC=2,
                     threshold_fdr=0.1) {

  # input should already have groups
  groups <- unique(sce$group)
  
  # calculate detection rate per gene
  cdr.dt <- data.table(
    gene = rownames(sce),
    detection_rate_A = DelayedArray::rowMeans(counts(sce[,sce$group==groups[1]])>0),
    detection_rate_B = DelayedArray::rowMeans(counts(sce[,sce$group==groups[2]])>0)
  ) %>% data.table::setnames(c("ens_id",sprintf("detection_rate_%s",groups[1]),sprintf("detection_rate_%s",groups[2])))


  # Filter genes by min detection rate per group
  cdr_A <- DelayedArray::rowMeans(counts(sce[,sce$group==groups[1]])>0) >= min_detection_rate_per_group
  cdr_B <- DelayedArray::rowMeans(counts(sce[,sce$group==groups[2]])>0) >= min_detection_rate_per_group
  sce <- sce[cdr_B | cdr_A,]


  # Convert SCE to DGEList
  sce_edger <- scran::convertTo(sce, type="edgeR")

  # Define design matrix (with intercept)
  cdr <- DelayedArray::colMeans(counts(sce)>0)
  design <- model.matrix(~ cdr + sce$batch + sce$group)

  # Estimate dispersions
  sce_edger  <- edgeR::estimateDisp(sce_edger, design)

  # Fit GLM
  fit <- edgeR::glmQLFit(sce_edger, design, coef=paste0('sce$group', groups[2]))

  # Test
  lrt <- edgeR::glmQLFTest(fit, coef=paste0('sce$group', groups[2]))

  # Construct output data.frame
  out <- edgeR::topTags(lrt, n=nrow(lrt))$table


  return(out)

}


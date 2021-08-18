
#'
#' @export
computePseudocells <- function(sce, assay="counts", order_by, kernel_size = 0.1,
                               num_cells = 200, density_thresh = 1) {

  exp_dense <- as.matrix(assay(sce,assay))
  pseudotime <- sce[[order_by]]

  #generate equally-spaced points along the trajectory:
  pseudo_locs <- seq(from=min(pseudotime), to=max(pseudotime),
                    length.out = num_cells)

  # remove pseudocells in low density areas
  dens <- density(pseudotime,bw=0.02)
  pseudo_dens <- approx(dens$x,dens$y,xout=pseudo_locs)
  pseudo_locs <- pseudo_locs[pseudo_dens$y > density_thresh]

  pseudo_list  <- lapply(1:length(pseudo_locs), function(i){
    pseudo_cell <- pseudo_locs[i]
    pseudo_dist <- pseudotime - pseudo_cell #length - samples number
    weights <- exp(-(pseudo_dist^2)/(kernel_size^2))
    norm_weights <- weights/sum(weights)
    out <- list(pseudo_mat=(exp_dense %*% norm_weights), pseudo_weights=weights)
    return(out)
  })

  pseudo_mat <- mapply(`[[`, pseudo_list, 1)
  pseudo_weights <- mapply(`[[`, pseudo_list, 2)

  return(list(pseudo_locs=pseudo_locs, pseudo_mat=pseudo_mat,
              pseudo_weights=pseudo_weights))
}



computePseudoEmbedding <- function(sce, weights, dimred="UMAP") {

  embedding <- as.matrix(reducedDim(sce,dimred))

  df <- t(apply(weights, 2, function(x) {
    pos <- (x%*% embedding)/sum(x)
    return(pos)
  }))

  return(df)
}

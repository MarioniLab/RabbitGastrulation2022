

#'
#'@export
plotPseudocells <- function(sce, pseudo_embed, dimred="UMAP",
                            colour_by="celltype",colours=celltype_colours,
                            point_size=1) {

  df <- reducedDim(sce,dimred)
  df <- as.data.frame(df)
  df <- cbind(df, sce[[colour_by]])
  colnames(df) <- c("x","y",colour_by)

  pseudo_embed <- as.data.frame(pseudo_embed)
  colnames(pseudo_embed) <- c("x","y")


  p <- ggplot2::ggplot(df, aes_string(x="x",y="y",colour=colour_by),
               size=point_size) +
    scale_color_manual(values = colours, name = "") +
    ggrastr::geom_point_rast(alpha=0.8,stroke=0,shape=16) +
    geom_point(data=pseudo_embed,aes(x=x,y=y),color="black",
               fill='red',stroke =0.5,size=point_size*1.5,shape=21) +
    theme_classic() +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    xlab(paste0(dimred,"_1")) + ylab(paste0(dimred,"_2")) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size=14),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          aspect.ratio = 1)


  return(p)

}





#' @importFrom stats aggregate
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggalluvial geom_alluvium geom_stratum geom_flow after_stat
#' @importFrom ggplot2 ggplot geom_text scale_x_discrete theme
#' @imporFrom ggplot2 scale_fill_manual theme_minimal aes element_blank
#' @importFrom ggplot2 element_text
#' @export
plotAnnotationAlluvium <- function(sce, old, new, ncell_thresh=10,
                                   repel_labels=TRUE, stratum_labels=FALSE,
                                   nudge=0.2) {

  df_frac <- aggregate(x = colnames(sce), by = list(sce[[old]],sce[[new]]),
                       FUN = length)
  colnames(df_frac) <- c("original", "predicted", "num_cells")
  df_frac[["original"]] <- as.factor(paste0("orig_",df_frac[["original"]]))
  df_frac[["predicted"]] <- as.factor(paste0("pred_", df_frac[["predicted"]]))


  # Filter entries according to the number of cells
  df_frac <- df_frac[df_frac$num_cells > ncell_thresh,]

  cols <- getCelltypeColours()
  names(cols) <- paste0("pred_",names(cols))


  p <- ggplot(df_frac,aes(axis1 = original, axis2 = predicted, y=num_cells, fill=predicted, label=predicted)) +
    geom_alluvium() +  geom_stratum() +
    scale_x_discrete(limits = c("Original annotation", "Predicted annotation"),
                     expand = c(.4, 0.4)) +
    geom_flow()

  if(repel_labels) {
    p <- p + ggrepel::geom_text_repel(
      aes(label = ifelse(substr(as.character(after_stat(stratum)), start = 1, stop = 5) == "orig_",
                         substring(as.character(after_stat(stratum)), 6) , NA)),
      stat = "stratum", size = 4, direction = "y", nudge_x = -nudge
    ) +

      ggrepel::geom_text_repel(
        aes(label = ifelse(substr(as.character(after_stat(stratum)), start = 1, stop = 5)=="pred_",
                           substring(as.character(after_stat(stratum)), 6), NA)),
        stat = "stratum", size = 4, direction = "y", nudge_x = nudge
      )
  }

  if(stratum_labels) {
    p <- p + geom_text(stat = "stratum",
                       aes(label = substring(as.character(after_stat(stratum)), 6)),
                       min.y=30)
  }

  p <- p + theme_minimal() +
    ylab("Number of cells") +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.background = element_blank(),
      text = element_text(size=20)
    ) +
    scale_fill_manual(values=cols) +
    theme(legend.position = "none")

  return(p)

}



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

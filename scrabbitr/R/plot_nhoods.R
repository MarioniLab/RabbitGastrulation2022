


#' @importFrom ggrastr rasterise
#' @importFrom miloR nhoodGraph
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom igraph vertex_attr
#' @importFrom ggraph geom_edge_link0 geom_node_point
#' @importFrom viridis scale_color_viridis
#' @importFrom ggplot2 theme theme_classic
plotNhoodMaxSim <- function(milo, df_maxNhood, colour_by="sim", legend_title="Max correlation" ) {

  nh_graph <- miloR::nhoodGraph(milo)
  
  # TODO: check df_maxNhood order is same as nh_graph
  V(nh_graph)$max_correlation <- df_maxNhood[[colour_by]]
  
  # data.frame / data.table implementations if order different
  #V(nh_graph)$max_correlation <- df_maxNhood[vertex_attr(nh_graph)$name, colour_by]
  #V(nh_graph)$max_correlation <- df_maxNhood[match(vertex_attr(nh_graph)$name, df_maxNhood[[id.col]]),]
  
  
  layout <- reducedDim(milo, "UMAP")[as.numeric(vertex_attr(nhoodGraph(milo))$name),]
  colnames(layout) <- c("x","y")

  p <- ggraph(simplify(nh_graph), layout = layout) +
    ggrastr::rasterise(geom_edge_link0(edge_colour = "grey66", edge_alpha=0.1),dpi=300) +
    ggrastr::rasterise(geom_node_point(aes(color = max_correlation), size = 3,alpha=0.8, shape=20),dpi=300) +
    theme_classic(base_size=14) +
    theme(axis.line = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), axis.title = element_blank(),
          aspect.ratio=1) +
    viridis::scale_color_viridis(option="cividis",name="Max correlation")

  return(p)

}



#' sim_values needs to be in order of milo nhoods
#' TODO: Split into smaller functions
plotNhoodSimGroups <- function(milo, sim_values, group_by="celltype",facet_by=NULL, type="ridge", orientation="vertical",colour_by=NULL,size=0.2,
                               subset=NULL,decreasing=FALSE,rel_min_height=0, show_rank=FALSE, rank_size=2, group_colours=scrabbitr::celltype_colours, xlabel="Cell type", ylabel="Correlation") {
  graph <- miloR::nhoodGraph(milo)

  df <- data.frame(nhood=vertex_attr(graph)$name,
                   max_sim=sim_values,
                   group=colData(milo)[as.numeric(vertex_attr(graph)$name), group_by])

  if(!is.null(facet_by)) {
    df$facet_by <- colData(milo)[as.numeric(vertex_attr(graph)$name), facet_by]
  }

  if(!is.null(colour_by)) {
    df$colour_by <- colData(milo)[as.numeric(vertex_attr(graph)$name), colour_by]
  }

  df_mean <- aggregate(df[,"max_sim"], list(df$group), mean,)
  df$group <- factor(df$group, levels=df_mean[order(df_mean$x,decreasing=decreasing),"Group.1"])

  df_mean$ranking[order(df_mean$x, decreasing=TRUE)] <- paste0("(",1:nrow(df_mean),"/", nrow(df_mean), ")")
  rownames(df_mean) <- df_mean$Group.1
  df$ranking <- df_mean[as.character(df$group), "ranking"]

  if(!is.null(subset)) {
    show_rank <- TRUE
    df <- df[df$group %in% subset,]
    df_mean <- df_mean[df_mean$Group.1 %in% subset,]
  }

  p <- ggplot(df, aes(x=max_sim,y=group,fill=group,color=colour_by)) + theme_bw()

  if(type=="boxplot") {
    p <- p + geom_boxplot(outlier.size=0.1)

  } else if(type=="violin") {
    p <- p + geom_violin()

  } else {
    p <- p + geom_density_ridges(size=size,rel_min_height = rel_min_height) +
      theme_ridges(grid = TRUE,center_axis_labels = TRUE)

    if(show_rank) {
      p <- p + geom_text(data=df_mean, aes(x=max(df$max_sim) + 0.05, y=Group.1, label=ranking, color=),
                         size = rank_size, inherit.aes = FALSE) +
        coord_cartesian(xlim = c(min(df$max_sim), max(df$max_sim)),clip = 'off')
    }
  }

  if(!is.null(facet_by)) p <- p + facet_grid(cols=vars(facet_by))

  if(orientation=="horizontal") {
    p <- p + coord_flip()
    if(!is.null(facet_by)) p <- p + facet_grid(rows=vars(facet_by))
  }

  p <- p +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_manual(values = group_colours[names(group_colours) %in% unique(df$group)], name = "") +
    guides(fill="none",color=guide_legend(title="",override.aes = list(fill = "white"))) +
    xlab(xlabel) + ylab(ylabel)

  return(p)

}






plotNhoodMappings <- function(r_milo, m_milo, df_sim, dimred="UMAP", colour_by,
                         r_graph=NULL, m_graph=NULL, r_umap= FALSE, m_umap=FALSE,
                         offset=c(10,0), reflect.X=FALSE, reflect.Y = FALSE,
                         rotate=NULL, line_alpha=0.05, edge_alpha=0.2,
                         colours=celltype_colours, legend_pos="none") {

  # Check if graph provided
  if(is.null(r_graph)) { r_graph <- nhoodGraph(r_milo) }
  if(is.null(m_graph)) { m_graph <- nhoodGraph(m_milo) }

  # Get nhood embedding positions
  r_nhoodPos <- getNhoodPositions(r_milo, r_graph)
  m_nhoodPos <- getNhoodPositions(m_milo, m_graph)

  # TODO: Fix this
  if(!is.null(rotate)) {
    rotate <- rotate*pi/180
    com <- colMeans(m_nhoodPos[,1:2])
    m_nhoodPos = sweep(m_nhoodPos,2,com, "-")
    m_nhoodPos[,1] = m_nhoodPos[,1]*cos(rotate) - m_nhoodPos[,2]*sin(rotate)
    m_nhoodPos[,2] = m_nhoodPos[,1]*sin(rotate) + m_nhoodPos[,2]*cos(rotate)
    m_nhoodPos = sweep(m_nhoodPos,2,com, "+")
  }
  if(reflect.X) { m_nhoodPos[,1] = -m_nhoodPos[,1]}
  if(reflect.Y) {  m_nhoodPos[,2] = -m_nhoodPos[,2]}


  # Add spacing to offset mouse and rabbit umaps
  m_nhoodPos[,1] <- m_nhoodPos[,1] + offset[1]
  m_nhoodPos[,2] <- m_nhoodPos[,2] + offset[2]

  r_df <- cbind(r_nhoodPos, vertex_attr(r_graph)[[colour_by]])
  colnames(r_df) <- c("x","y","obs")

  m_df <- cbind(m_nhoodPos, vertex_attr(m_graph)[[colour_by]])
  colnames(m_df) <- c("x","y","obs")

  # Plot nhoods and nhood graphs
  df_plot <- rbind(r_df, m_df)
  p <- ggplot(df_plot,aes(x=x,y=y)) +
    geom_edge_link0(data=get_edges(format="short")(create_layout(simplify(r_graph),layout=r_nhoodPos)),
                    edge_colour = "grey66", edge_alpha=edge_alpha) +
    geom_edge_link0(data=get_edges(format="short")(create_layout(simplify(m_graph),layout=m_nhoodPos)),
                    edge_colour = "grey66", edge_alpha=edge_alpha) +
    geom_point(aes(fill=obs),stroke=0,shape=21) +
    scale_fill_manual(values = celltype_colours[names(celltype_colours) %in% unique(df_plot$obs)], name = "")


  # Filter sim data frame to nhoods within graphs provided
  r_nhoodIDs <- as.numeric(vertex_attr(r_graph)$name)
  m_nhoodIDs <- as.numeric(vertex_attr(m_graph)$name)

  sim_filt <- df_sim[(df_sim[,1] %in%  r_nhoodIDs) &
                       (df_sim[,2] %in% m_nhoodIDs),]


  # Give each line a unique name
  sim_filt$alignment <- paste0("align_",sim_filt$r_nhood,"_",sim_filt$m_nhood)

  # Link rabbit and mouse nhood positions by alignment name
  r_lines <- sim_filt
  r_lines[,c("x","y")] <- r_nhoodPos[as.character(r_lines[,1]),]

  m_lines <- sim_filt
  m_lines[, c("x","y")] <- m_nhoodPos[as.character(m_lines[,2]),]
  df_lines <- rbind(r_lines, m_lines)


  # Add similarities lines
  p <- p + geom_line(data=df_lines,aes(x=x,y=y, group=alignment,colour=sim),alpha=line_alpha) +
    scale_color_distiller(palette = "Spectral",limits=c(min(df_sim$sim),max(df_sim$sim)))


  # Add theme
  p <- p +  theme_classic(base_size=14) +
    theme(axis.line = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), axis.title = element_blank(),
          legend.position=legend_pos,aspect.ratio=1)


  return(p)

}



plotTrajMappings <- function(r_milo, m_milo, df_sim, group_by, groups,
                             dimred = "UMAP", ...) {

  r_traj <- subsetMiloGroups(r_milo, group_by, groups)
  m_traj <- subsetMiloGroups(m_milo, group_by, groups)

  p <- plotNhoodMappings(r_traj, m_traj, df_sim, dimred, group_by, ...)

  return(p)

}



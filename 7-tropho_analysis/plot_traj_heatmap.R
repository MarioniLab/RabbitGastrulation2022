
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(SingleCellExperiment)
library(RColorBrewer)
library(data.table)

setwd("G:/My Drive/Postgrad/PhD/Projects/RabbitGastrulation2022/")

# Load data
r_syn <- readRDS("data-in/trophoblast/r_syn.rds")

# Rename cell types
r_syn$celltype[r_syn$celltype == "Syncytiotrophoblast progenitors"] <- "SCT progenitors"
r_syn$celltype[r_syn$celltype == "Syncytiotrophoblast"] <- "Early SCT"


plot_genes <- c("CDX2", "TEAD4", "HAND1","DLX5","DLX6","ZNF740","PATZ1",
                "TFEB", "GCM1","CEBPA","ELF3","TFAP2C", "MITF","TP63")

tropho_celltypes <- c("Trophoblast", "Cytotropoblast", 
                      "SCT progenitors", 
                      "Early SCT")

hm_width <- unit(5.5, "cm")

n <- 49

moving_average <- function(x, n = 49) {
  ret  <- cumsum(x)
  ret[(n+1):(length(ret))] = ret[(n+1):(length(ret))] - ret[1:(length(ret)-n)]
  return(ret[(n) :length(ret)] / n)
}

minmax <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}

values <- as.matrix(logcounts(r_syn)[plot_genes,order(r_syn$dpt_pseudotime)])
colnames(values) <- NULL

pad_size <- (n-1)/2
ma_window <- (n/2):(ncol(r_syn) - (n/2))
padding <- matrix(0, ncol=pad_size, nrow = nrow(values))
pad_values <- cbind(padding,values,padding)

ma_values <- t(apply(values,1, moving_average, n=n))
norm_values <- t(apply(ma_values, 1, minmax))

pseudo_col <- colorRamp2(seq(0, 1, by = 0.25), magma(5))
rna_annot <- HeatmapAnnotation(celltype = r_syn$celltype[order(r_syn$dpt_pseudotime)][ma_window],
                               stage = r_syn$stage[order(r_syn$dpt_pseudotime)][ma_window],
                               pseudotime = r_syn$dpt_pseudotime[order(r_syn$dpt_pseudotime)][ma_window], 
                               col = list(celltype = scrabbitr::getCelltypeColours(),
                                          stage = scrabbitr::getStageColours(),
                                          pseudotime = pseudo_col),
                               annotation_name_side = "left")



colours <- scico::scico(3,palette = "vik")
rna_col <- colorRamp2(c(0,1), colours[2:3])

p1 <- Heatmap(norm_values, use_raster = TRUE,
              cluster_columns = FALSE,
              cluster_rows =FALSE,
              bottom_annotation = rna_annot,
              col = rna_col,
              row_names_side = "left",
              width = hm_width,
              column_title ="RNA-seq",
              name = "Norm exp\n(min-max)"
              )


atac_raw <- fread(file.path("data-in/trophoblast/atac_data", 'trophoblast_atac_motifs.csv'), 
                  data.table = FALSE)
rownames(atac_raw) <- atac_raw[,1]
atac_raw <- atac_raw[,-1]
atac_raw <- atac_raw[plot_genes,]


atac_meta <- fread(file.path("data-in/trophoblast/atac_data", 'trophoblast_atac_meta.csv'), 
                   data.table = FALSE)

atac_meta <- atac_meta[order(atac_meta$Pseudotime),]
atac_meta$celltype_tropho[atac_meta$celltype_tropho == "Syncytiotrophoblast Progenitors"] <- "SCT progenitors"
atac_meta$celltype_tropho[atac_meta$celltype_tropho == "Syncytiotrophoblast"] <- "Early SCT"
atac_meta$stage[atac_meta$stage == "GD9_ExE"] <- "GD9"

atac_raw <- atac_raw[,atac_meta$cell]
atac_raw <- atac_raw[,atac_meta$cell[!is.na(atac_meta$Pseudotime)]]
atac_meta <- atac_meta[!is.na(atac_meta$Pseudotime),]
colnames(atac_raw) <- NULL

ma_window <- (n/2):(ncol(atac_raw) - (n/2))
padding <- matrix(0, ncol=pad_size, nrow = nrow(atac_raw))
pad_values <- cbind(padding,atac_raw,padding)

ma_values <- t(apply(atac_raw,1, moving_average, n = n))
norm_values <- t(apply(ma_values, 1, scale))

pseudo_col <- colorRamp2(seq(0, 100, by = 25), magma(5))
atac_annot <- HeatmapAnnotation(celltype = atac_meta$celltype_tropho[ma_window],
                                stage = atac_meta$stage[ma_window],
                                pseudotime = atac_meta$Pseudotime[ma_window], 
                                col = list(celltype = scrabbitr::getCelltypeColours(),
                                           stage = scrabbitr::getStageColours(),
                                           pseudotime = pseudo_col),
                                show_annotation_name  = FALSE)

atac_col <- colorRamp2(c(min(norm_values), 0, max(norm_values)), 
                       colours)

p2 <- Heatmap(as.matrix(norm_values), use_raster = TRUE,
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              bottom_annotation = atac_annot,
              column_split = ifelse(atac_meta$stage[ma_window] == "GD7","1","2"),
              cluster_column_slices = FALSE,
              col = atac_col,
              row_names_side = "left",
              width = hm_width,
              column_title = "ATAC-seq",
              name = "Norm motif\naccessibility\n(z-score)")



p <- p1 + p2
draw(p, ht_gap = unit(0.5, "cm"))



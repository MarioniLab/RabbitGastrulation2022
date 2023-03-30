
library(SingleR)

# Load data
r_sce <- readRDS("data-in/r_sce.rds")
m_sce <- readRDS("data-in/m_sce.rds")
rm_orthologs <- read.table("data-in/rm_orthologs.tsv",sep="\t")


# Train SingleR model
m_model <- trainSingleR(ref = m_sce[rm_orthologs$mm,], labels=m_sce$celltype,
                       de.method="wilcox", de.n=50, aggr.ref = T)

# Save model
saveRDS(m_model, "data-out/m_singler_model.rds")

# Load model
# m_model <- readRDS("data-out/m_singler_model.rds")


# Ensure common features
r_one2one <- r_sce[rm_orthologs$oc,]
rownames(r_one2one) <- sfs$mm


# Classify rabbit cells
r_pred <- classifySingleR(test=r_one2one, trained = m_model) 

# Export results
saveRDS(r_pred, "data-out/r_singler_out.rds")



# Normalise correlation scores (used SingleR::plotScoreHeatmap)
mmax <- rowMaxs(sr_model$scores, na.rm = TRUE)
mmin <- rowMins(sr_model$scores, na.rm = TRUE)
scores <- (sr_model$scores-mmin)/pmax(mmax-mmin, 1e-8)
scores <- scores^3


# Plot score heatmap (Supplementary Figure 4B)
ha <- HeatmapAnnotation(singler = r_data$singler,
                       celltype = r_data$celltype,
                       col = list(celltype = cmap[r_data$celltype],
                                  singler = cmap[r_data$singler]), 
                       show_legend = FALSE)

#pdf(file="../../plots/revisions/singler_scores_ordered.pdf", width = 10 , height = 14)
ht <- Heatmap(t(scores),
              name = "Score",
              cluster_rows = F,
              cluster_columns = FALSE,
              show_column_names = FALSE,
              column_order = order(r_data$singler),
              top_annotation = ha,
              col = circlize::colorRamp2(c(0,0.5,1),
                                         viridis::viridis(3,option="D")),
              use_raster = TRUE,
              raster_quality = 5)

draw(ht)
#dev.off()





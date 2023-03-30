
library(SingleCellExperiment)
library(scran)
library(irlba)
library(BiocParallel)
library(irlba)
library(batchelor)

ncores <- 4

r_sce <- readRDS("data-in/r_data_preprocessed.rds")

# Get HVGs
gene_var <- modelGeneVar(r_sce)
hvgs <- getTopHVGs(gene_var, n=3000)

# Compute n_cells in each sample
meta <- colData(r_sce)
df <- meta[, c("stage", "sample")]
df$ncells = sapply(df$sample, function(x) sum(meta$sample == x))
df$stage = factor(df$stage, levels = rev(c("GD9","GD8","GD7")))
df$sample <- as.character(df$sample)
df$stage = as.character(df$stage)

# Order by number of cells
df = df[order(df$stage, df$ncells, decreasing = TRUE),]
merge_order <- rev(lapply(split(df,df$stage), function(x) list(unique(x$sample))))


# Combine GD9 samples by anatomical location
gd9_head <- list("12","16")
gd9_trunk <- list("13","17")
gd9_tail <- list("14","18")
gd9_ep <- list(list("23","24","25"),list("20","21"))
gd9_ys <- list("15","19","22","26")

merge_order[[1]] <- list(gd9_ep,gd9_ys,gd9_trunk,gd9_tail,gd9_head)
names(merge_order)<- NULL


# Run fastMNN
out <- fastMNN(r_sce,batch=r_sce$sample,k=20,subset.row = hvgs,
               merge.order = merge_order,BPPARAM = MulticoreParam(ncores))

# Export results
saveRDS(out,"data-out/fastMNN_out.rds")
write.table(reducedDim(out,"corrected"),"data-out/corrected_pcs.tsv",sep="\t",quote=F)


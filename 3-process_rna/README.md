# scRNA-seq processing

## Methods

### Cell calling

Cell calling and doublet removal was adapted from (Pijuan-Sala et al. 2019). Briefly, cell barcodes that were associated with real cell transcriptomes were identified using emptyDrops (Lun et al. 2019), which assesses cells with RNA content distinct from ambient background RNA, the latter determined from barcodes associated with fewer than 100 unique molecular identifiers (UMIs). Cells with P \< 0.01 (Benjamini--Hochberg-corrected) and at least 3500 UMIs and 900 unique genes were considered for further analysis.

Additionally, cells with mitochondrial gene-expression fractions greater than 24.03% were excluded. Rabbit mitochondrial fractions were suspected to be higher than mouse due to species differences or genome annotation. The thresholds were determined by considering a median-centered median absolute deviation (MAD)-variance normal distribution; cells with mitochondrial read fraction outside of the upper end of this distribution were excluded (where outside corresponds to P \< 0.05; Benjamini--Hochberg-corrected).

### Doublet detection

Doublets were scored as previously described in Pijuan-Sala et al. 2019. First, a doublet score was computed for each cell by applying the `doubletCells` function (scran R package) to each sample separately. This function returns the density of simulated doublets around each cell, normalized by the density of observed cell libraries. High scores indicate high doublet probability. We next identified clusters of cells in each sample by computing the first 50 principal components (PCs) across all genes, building a shared nearest-neighbour graph (10 nearest neighbours; `buildSNNGraph` function; scran R package), and applying the Louvain clustering algorithm (`cluster_louvain` function; igraph R package; default parameters) to it. Only HVGs (calculated separately for each sample) were used for the clustering. This procedure was repeated in each identified cluster to break the data into smaller clusters, ensuring that small regions of high doublet density were not clustered with large numbers of singlets. For each cluster, the median doublet score was considered as a summary of the scores of its cells, as clusters with a high median score were likely to contain mostly doublets. Doublet calls were made in each sample by considering a null distribution for the scores using a median-centred MAD-variance normal distribution, separately for each sample. The MAD estimate was calculated only on values above the median to avoid the effects of zero-truncation, as doublet scores cannot be less than zero. All cells in clusters with a median score at the extreme upper end of this distribution (Benjamini--Hochberg-corrected P \< 0.1) were labelled as doublets. A final clustering step was performed across all samples together to identify cells that shared transcriptional profiles with called doublets, but escaped identification in their own samples. Clusters were defined using the same procedure as was applied to each sample, with the exceptions that sub-clustering was not performed, and batch-corrected principal components were used (see 'Batch correction', above). To identify clusters that contained more doublets than expected, we considered for each cluster the fraction of cell libraries that were called as doublets in their own samples. We modeled a null distribution for this fraction using a median-centered, MAD-estimated variance normal distribution as described for the median doublet score in each sample, above, and called doublets from the distribution as in each sample, above.

### Normalisation

Transcriptome size factors were calculated as previously described (Pijuan-Sala et al. 2019) using `computeSumFactors` from the scran R package (version 1.18.7). Cells were pre-clustered with the `quickCluster` function using the parameter `method=igraph` (using the scran R package), and minimum and maximum cluster sizes of 100 and 3,000 cells, respectively. Raw counts for each cell were divided by their size factors, and the resulting normalized counts were used for further processing.

## References

See `references.bib`.

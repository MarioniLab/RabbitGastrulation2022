# Compare species

This directory contains code to perform SAMap integration and run the neighbourhood comparison pipeline to compare cell states across rabbit and mouse scRNA-seq data.

`run_samap.ipynb` - Performs SAMap integration to jointly embed the rabbit and mouse datasets.

`compare_nhoods.ipynb` - Runs through the neighbourhood comparison pipeline.

## Requirements

The code deposited here assumes that the rabbit scRNA-seq data has been processed and annotated. Please see the other directories in this repository to reproduce these previous processing steps. Fully-processed forms of the rabbit scRNA-seq atlas can be obtained from [here](https://content.cruk.cam.ac.uk/jmlab/RabbitGastrulation2022/data/RNA/).

These notebooks compare the rabbit atlas with the extended mouse atlas (Imaz Rosshandler et al. 2023) available from <https://marionilab.github.io/ExtendedMouseAtlas/#data>.

The neighbourhood comparison pipeline and associated plotting functions are implemented in the [scrabbitr](https://github.com/dkeitley/scrabbitr) package. This can be installed from the Github repository: <https://github.com/dkeitley/scrabbitr>. `scrabbitr` is heavily dependent on the [miloR](https://github.com/MarioniLab/miloR) package.

SAMap can be installed from the Github repository: <https://github.com/atarashansky/SAMap>. Note that in this work, we used version 0.1.6. Please refer to the SAMap Github repository for more up to date vignettes/documentation.

## Additional Methods

### Neighbourhood pipeline

To compute the correlation between neighbourhoods, we used the intersection of the top 2000 highly variable genes that are one-to-one orthologs across the two datasets. The intersection of genes was chosen to avoid confounding differences in expression with technical variation resulting from the lower quality rabbit genome annotation. It ensured that genes selected for comparisons were expressed and highly variable across both datasets. We also experimented with the number of highly variable genes but found that our results changed very little above 2000 HVGs.

In Figures 4 and Extended Data Figure 7, we visualised neighbourhood similarities across the rabbit and mouse UMAP embeddings. For each rabbit neighbourhood we extracted its maximum similarity score with any mouse neighbourhood (and vice versa). We then plotted each neighbourhood according to the UMAP position of its index cell, and coloured each point by the neighbourhood's maximum correlation value. In Figures 4D, the maximum correlation values are aggregated according to the cell type annotation of each index cell, to obtain a distribution of similarity scores across each cell type. Finally, for the trajectory comparisons, subsets of neighbourhoods were extracted from both the rabbit and mouse datasets whose index cell was annotated as one of a specified set of cell types. These neighbourhoods were then plotted in the same way as Figure 4C and repositioned or reflected to facilitate a clear visual comparison. Lines were drawn between maximally correlated neighbourhood pairs, computed in both the rabbit-mouse and mouse-rabbit directions. The line colour indicates the strength of correlation for each mapping.

### Table of Contents

1. [Paper abstract](#Abstract)

2. [Data availability](#Data availability)

3. [Explore the data](#Explore the data)

4. [Code availability](#Code availability)

5. [Contact](#Support or Contact)

6. [Other links](#Other links)

   

## Rabbit Development as a Model for Single Cell Comparative Genomics



**Mai-Linh Ton<sup>1,2,\*</sup> Daniel Keitley<sup>3,*</sup>, Bart Theeuwes<sup>2</sup>, Carolina Guibentif<sup>4</sup>, Jonas Ahnfelt-Rønne<sup>5</sup>, Ivan Imaz-Rosshandler<sup>6</sup>, Èlia Benito-Gutiérrez<sup>3</sup>, John Marioni<sup>7,8,9</sup>, Berthold Göttgens<sup>1,2</sup>** 



1. *Department of Haematology, University of Cambridge, Cambridge, UK*

2. *Wellcome-Medical Research Council Cambridge Stem Cell Institute, University of Cambridge, Cambridge, UK*

3. *Department of Zoology, University of Cambridge, Cambridge, UK*

4. *Department of Microbiology and Immunology, University of Gothenburg, Gothenburg, Sweden*

5. *Department of Pathology & Imaging, Novo Nordisk, Måløv, Denmark*

6. *Medical Research Council Laboratory of Molecular Biology, Cambridge, UK*

7. *Wellcome Sanger Institute, Wellcome Genome Campus, Cambridge, UK*

8. *European Molecular Biology Laboratory European Bioinformatics Institute, Cambridge, UK*

9. Cancer Research UK Cambridge Institute, University of Cambridge Cambridge, UK

   

### Abstract

Biomedical research relies heavily on the use of model organisms to gain insight into human health and development.  Traditionally, the mouse has been the favored vertebrate model, due to its experimental and genetic tractability. Non-rodent embryological studies however highlight that many aspects of early mouse development, including the egg-cylinder topology of the embryo and its method of implantation, diverge from other mammals, thus complicating inferences about human development. In this study, we constructed a morphological and molecular atlas of rabbit development, which like the human embryo, develops as a flat-bilaminar disc. We report transcriptional and chromatin accessibility profiles of almost 200,000 single cells and high-resolution histology sections from embryos spanning gastrulation, implantation, amniogenesis, and early organogenesis. Using a novel computational pipeline, we compare the transcriptional landscape of rabbit and mouse at the scale of the entire organism, revealing that extra-embryonic tissues, as well as gut and PGC cell types, are highly divergent between species. Focusing on these extra-embryonic tissues, which are highly accessible in the rabbit, we characterize the gene regulation underlying trophoblast differentiation and identify novel signaling interactions involving the yolk sac mesothelium during hematopoiesis. Finally, we demonstrate how the combination of both rabbit and mouse atlases can be leveraged to extract new biological insights from sparse macaque and human data. The datasets and analysis pipelines reported here will expedite the development of models of mammalian gastrulation and set a framework for a broader cross-species approach to decipher early mammalian development.



### Data availability

Various forms of the transcriptomics data are available [here](https://content.cruk.cam.ac.uk/jmlab/RabbitGastrulation2022/) for loading into R and python. 

| File name                                                    | Description                                                  |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| `r_data.h5ad`                                                | AnnData object for processing with [scanpy](https://scanpy.readthedocs.io/en/stable/index.html). |
| `r_sce.rds`                                                  | Contains a `SingleCellExperiment` object for processing in R. |
| `r_counts.mtx`                                               | Counts matrix in MatrixMarket format                         |
| `r_logcounts.mtx`                                            | Normalised logcounts in MatrixMarket format                  |
| `r_meta.tsv`                                                 | Per-cell observations (e.g. sample, cell type annotation, stage) |
| `r_genes.tsv`                                                | Ensembl codes and gene names for the scRNA-seq features.     |
| `r_sizefactors.tsv`                                          |                                                              |
| `r_corrected_pcs.tsv`                                        | The top 50 principal components resulting from fastMNN batch correction. |
| `r_reducedDims.rds` `r_umap.tsv` `r_tsne.tsv` `r_fa.tsv` `r_umap3d.tsv` | List of DataFrames containing coordinates for reduced dimensionality representations (e.g. UMAP, TSNE, PCA). |
| `oc_mm_orthologs.tsv`                                        | DataFrame of one-to-one orthologs between the rabbit (*Oryctolagus cuniculus*) and mouse (*Mus musculus*) obtained through Ensembl. |

Raw scRNA-seq files, ATAC-seq data and histology images will be made available in due course. For details of the other, externally generated datasets used in our analysis, see the methods section of the paper. 



### Explore the data

#### Shiny app / Vitessce

You can interactively explore our single-cell transcriptomic dataset on the the web via our shiny app accessible at [https://crukci.shinyapps.io/scrabbit-shiny/](https://crukci.shinyapps.io/scrabbit-shiny/). 

#### cellxgene

The data can also be explored locally using [cellxgene](https://github.com/chanzuckerberg/cellxgene). After following the cellxgene installation instructions, the rabbit data can be loaded via the anndata `.h5ad` file stored in the link above. 

```
cellxgene launch https://content.cruk.cam.ac.uk/jmlab/RabbitGastrulation2022/r_data.h5ad
```



### Code availability

The code used to process and analyse and the rabbit transcriptomics data is available through the [RabbitGastrulation2022](https://github.com/dkeitley/RabbitGastrulation2022) github repository.

The code is organised into R and python Jupyter notebooks. Many of the functions used are packaged into the [scrabbitr](https://github.com/dkeitley/scrabbitr) and [scrabbitpy](https://github.com/dkeitley/scrabbitpy) packages. See the individual github repositories for installation instructions. 



### Support or Contact

General queries can be directed to [Bertie Göttgens](bg200@cam.ac.uk) , [John Marioni](mailto:marioni@ebi.ac.uk) or [Èlia Benito-Gutiérrez](mailto:eb647@cam.ac.uk). For issues relating to the data, code or shiny app, you can file an issue on the most relevant github repository or email Daniel Keitley at [dk562@cam.ac.uk](mailto:dk562@cam.ac.uk). 



### Other links

[Göttgens lab website](https://www.stemcells.cam.ac.uk/people/pi/gottgens)

[Marioni lab website](https://www.ebi.ac.uk/research-beta/marioni/)

[Benito-Gutiérrez lab website](https://www.zoo.cam.ac.uk/research/cell-and-developmental-biology/benito-gutierrez)


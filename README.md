# RabbitGastrulation2022

This repository contains functions and scripts used as part of the Ton et al. 2022 paper on rabbit gastrulation. 



> **Rabbit Development as a Model for Single Cell Comparative Genomics**

>
> Mai-Linh Ton<sup>1,2,\*</sup> Daniel Keitley<sup>3,*</sup>, Bart Theeuwes<sup>2</sup>, Carolina Guibentif<sup>4</sup>, Jonas Ahnfelt-Rønne<sup>5</sup>, Thomas Kjærgaard Andreassen<sup>5</sup>, Fernando J. Calero-Nieto<sup>2</sup>, Ivan Imaz-Rosshandler<sup>6</sup>, Blanca Pijuan-Sala<sup>7</sup>, Jennifer Nichols<sup>2</sup>, Èlia Benito-Gutiérrez<sup>3</sup>, John Marioni<sup>7,8,9</sup>, Berthold Göttgens<sup>1,2</sup>

> **Abstract:**
Biomedical research relies heavily on the use of model organisms to gain insight into human health and development.  Traditionally, the mouse has been the favored vertebrate model, due to its experimental and genetic tractability. Non-rodent embryological studies however highlight that many aspects of early mouse development, including the egg-cylinder topology of the embryo and its method of implantation, diverge from other mammals, thus complicating inferences about human development. In this study, we constructed a morphological and molecular atlas of rabbit development, which like the human embryo, develops as a flat-bilaminar disc. We report transcriptional and chromatin accessibility profiles of almost 200,000 single cells and high-resolution histology sections from embryos spanning gastrulation, implantation, amniogenesis, and early organogenesis. Using a novel computational pipeline, we compare the transcriptional landscape of rabbit and mouse at the scale of the entire organism, revealing that extra-embryonic tissues, as well as gut and PGC cell types, are highly divergent between species. Focusing on these extra-embryonic tissues, which are highly accessible in the rabbit, we characterize the gene regulation underlying trophoblast differentiation and identify novel signaling interactions involving the yolk sac mesothelium during hematopoiesis. Finally, we demonstrate how the combination of both rabbit and mouse atlases can be leveraged to extract new biological insights from sparse macaque and human data. The datasets and analysis pipelines reported here will expedite the development of models of mammalian gastrulation and set a framework for a broader cross-species approach to decipher early mammalian development.



## Requirements

This project relies on both R and python code packaged into [scrabbitr](https://github.com/dkeitley/scrabbitr) and [scrabbitpy](https://github.com/dkeitley/scrabbitpy) respectively. We recommend installing these packages within R and python  [conda environments](https://docs.conda.io/projects/conda/en/latest/index.html). The environments originally used in our analysis can be found in the [envs](envs/) directory.  

##### Installing scrabbitr

```R
# Inside R session
devtools::install_github("https://github.com/dkeitley/scrabbitr")
```

##### Installing scrabbit

```bash
git clone https://github.com/dkeitley/scrabbitpy.git
cd scrabbitpy
pip install .
```


## Usage / data availability
Each Jupyter notebook loads and exports data relevant to each step of analysis. 

Raw sequencing files are stored in ArrayExpress with accession numbers (scRNA-seq) **E-MTAB-11836** and (scATAC-seq) **E-MTAB-11804**.

Processed data can be obtained [here](https://content.cruk.cam.ac.uk/jmlab/RabbitGastrulation2022/). 

To explore these datasets interactively, see info [here](https://marionilab.github.io/RabbitGastrulation2022/#explore). 




## Support or Contact
General queries about the project can be directed to [Bertie Göttgens](bg200@cam.ac.uk) , [John Marioni](mailto:marioni@ebi.ac.uk) or [Èlia Benito-Gutiérrez](mailto:eb647@cam.ac.uk). For issues relating to the data, code or shiny app, you can file an issue on the most relevant github repository or email Daniel Keitley at [dk562@cam.ac.uk](mailto:dk562@cam.ac.uk). 




## Links
[Project webpage](https://marionilab.github.io/RabbitGastrulation2022/)

[Shiny app](https://crukci.shinyapps.io/scrabbit-shiny/)

[scrabbitr](https://github.com/dkeitley/scrabbitr)

[scrabbitpy](https://github.com/dkeitley/scrabbitpy)

[scrabbit-shiny](https://github.com/dkeitley/scrabbit-shiny)


[Göttgens lab website](https://www.stemcells.cam.ac.uk/people/pi/gottgens)

[Marioni lab website](https://www.ebi.ac.uk/research-beta/marioni/)

[Benito-Gutiérrez lab website](https://www.zoo.cam.ac.uk/research/cell-and-developmental-biology/benito-gutierrez)







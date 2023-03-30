# Annotate cell types

This directory contains code used to perform an automated and manual annotation of cell types in the rabbit scRNA-seq data.

`run_singler.R` - R script used to run the SingleR automated annotation.

`assign_cell_types.ipynb` - Python notebook where clusters are assigned cell type identities.

## Requirements

Data from the extended mouse atlas (Imaz Rosshandler et al. 2023) was used to train the automated annotation model. This can be accessed from <https://marionilab.github.io/ExtendedMouseAtlas/#data>.

`assign_cell_types.ipynb` uses utility functions defined in `scrabbitpy`. This can be installed from the Github repository: <https://github.com/dkeitley/scrabbitpy>

The scripts used here assume that the rabbit atlas has been processed, clustered and visualised using UMAP dimensionality reduction. To replicate these previous steps, please see the other directories in this repository.

## Further details

### Relation to the ExtendedMouseAtlas cell types

47/67 cell types defined in the rabbit were annotated consistently with the mouse. Based on evidence from marker gene expression, automated cell type annotation and SAMap integration, these cell types presented a 1-1 mapping with those in the mouse (Extended Data Figure 4A,C). Cell types were annotated differently between the two datasets if they were i) not confidently identifiable in the rabbit or mouse, ii) if a different terminology was preferred or iii) if a less-specific annotation was more suitable.

#### Less specific annotations

5/67 cell types fell into this last category - particularly neural cell types. Forebrain, midbrain and hindbrain populations were detectable in the rabbit atlas, based on signatures of expression (Forebrain: OTX2+, FEZF1/2+, SIX6+; Midbrain: OTX2+, FEZF2-, HOXA2-; Hinbrain: HOXA2+, HOXA7-, HOXB9-, OTX2-). However, the mouse forebrain was subclassified into *'Ventral forebrain progenitors'*, *'Late dorsal forebrain progenitors'* and *'Early dorsal forebrain progenitors'*, which are not as clearly defined in our rabbit data. Similarly the mouse atlas distinguishes *'Midbrain progenitors'* from *'Dorsal midbrain neurons'* and specifies a *'Midbrain/Hindbrain boundary'*. These separations were inconspicuous in the rabbit atlas and so a lower resolution classification was used. The same was true of the *'Hindbrain'* (subclassified into *'Hindbrain neural progenitors', 'Ventral hindbrain progenitors'* and *'Dorsal hindbrain progenitors'* in the mouse) and *'Spinal cord'* (annotated as *'Dorsal spinal cord progenitors'* and *'Spinal cord progenitors'* in the mouse). Mouse *'Branchial arch neural crest'* and *'Frontonasal mesenchyme'* were also grouped into a *'Cranial neural crest'* cell type in the rabbit. The automated label transfer of these cell types did not form distinct clusters in the rabbit, although they were proximally located in the UMAP embedding (Extended Data Figure 4A) and aligned closely with cells from anterior GD9 samples (Figure 2B).

#### Differences in terminology

5/67 cell types were annotated inconsistently due to different choices of terminology. For example, the rabbit cluster of *'Differentiating neurons'* exhibited strong similarity with *'Dorsal midbrain neurons'* of the mouse (Extended Data Figure 4A,C), however, there was no strong evidence indicating the spatial specificity of these cells in the rabbit, suggesting that a less-specific annotation was more suitable. These cells differentially express markers of early neural differentiation (NEUROD4 and ASCL1). The spatial specificity of the rabbit *'Floor plate'* (SHH+) was also indiscernible and so was labelled more generally than the *'Hindbrain floor plate'* annotation in the mouse. The neighbourhood comparisons and SAMap integration suggest that the rabbit rA3 cluster is closely related with the mouse *'Non-neural ectoderm 3'* cell type (Extended Data Figure Figure S4A, C, S8B). However these cells express markers of the amniotic ectoderm, such a GABRP and VTCN1 (Rostovskaya et al. 2022), and so this cluster was annotated as *'Amnion 3'* in the rabbit. *'Trophoblast'* and *'Hypoblast'* cell types were also renamed from *'ExE ectoderm'* and *'ExE endoderm'* respectively to reflect the uncertain relationship between extra-embryonic tissues and match terminology more commonly used in non-rodent embryology.

#### Cell type annotations specific to the rabbit atlas

Finally, 10/67 cell types, were not present in the mouse atlas and were annotated de-novo in the rabbit. This mostly includes extra-embryonic ectoderm cell types (*Amnion 1, Amnion 2, Cytotrophoblast, SCT progenitors, Early SCT*). Genes used to annotate these clusters are shown in Extended Data Figure Figure S5D. Concretely, the amnion 1 cluster was annotated based on the expression of markers associated with the first-wave of amniogenesis (Rostovskaya et al. 2022). This includes VTCN1, CDX2 (Rostovskaya et al. 2022) and TFAP2C (Yang et al. 2021), as well as aquaporin and prostaglandin genes, additionally implicated in amniotic fluid homeostasis (Martinez and Damiano 2017). Other amnion markers, KRT7 and WNT6 (Chuva de Sousa Lopes et al. 2022) were used to annotate *Amnion 2*. The third amnion cluster, *Amnion 3* specifically exhibited strong similarity with surface-ectoderm derived amniotic ectoderm and the second-wave of amniongenesis, expressing GABRB, ISL1 and WNT6 (Rostovskaya et al. 2022). These genes have been reported in the mouse, macaque and human amniotic ectoderm (Roost et al. 2015; Chuva de Sousa Lopes et al. 2022; Yang et al. 2021; Rostovskaya et al. 2022). Markers of the trophoectoderm, such as HAND1, CDX2 (Strumpf et al. 2005) and TEAD4 (Yagi et al. 2007) were used to annotate the *Trophoblast* which also primarily originates from the earliest GD7 timepoint (Figure 2B). The cytotrophoblast, consisting mostly of GD8 cells, was annotated based on a similar transcriptional signature (TEAD4+, HAND1+, WNT2,) as well as the expression of FGFR2 (D. Baczyk et al. 2006) and CYP19A1, which is upregulated during cytotrophoblast-syncytiotrophoblast differentiation (Kwak et al. 2019). *Synctiotrophoblast progenitors* were identified through the expression of GCM1 and CEBPA, key regulators of synctial fusion (Yu et al. 2002; D. Baczyk et al. 2009). Finally the *early synctiotrophoblast* population was annotated based on the expression of TFAP2A/C and TP63, which are similarly detected in human and mouse syncytiotrophoblast (Kuckenberg, Kubaczka, and Schorle 2012; Knofler et al. 2019). In addition to the extra-embryonic ectoderm clusters, the roof plate was identifiable as a distinct neural cluster that differentially expressed LMX1A (Millonig, Millen, and Hatten 2000). A definitive endoderm annotation was also absent from the mouse atlas and so this was annotated based on the expression of GSC, EOMES and FOXA2 (Lewis and Tam 2006). While the definitive endoderm was straightforward to detect, the relationship between more differenitated rabbit and mouse endoderm cell types was less clear, specifically cell types of the gut. The automated label transfer model annotated a trajectory of FOXA2+ cells as the *'Midgut'* (Extended Data Figure 4A). However, while the mouse defines foregut, midgut and hindgut cell types, these subdivisions are not identifiable in the rabbit and are potentially associated with later developmental stages (Extended Data Figure 8A). Hence the GD9 subset of predicted gut cells (DKK1+, PLAT+) were annotated as the *'Gut tube'*, (consistently with the mouse), while those predominately at the earlier GD8 timepoint was labelled as *'Gut endoderm'*. The gut tube annotations of the rabbit and mouse overlap in the SAMap embedding (Extended Data Figure S4C). *'Visceral YS endoderm'* and *'Visceral YS endoderm 2'* clusters were also given unique labels as their relationship with the mouse was not clear. These cells express markers of the extra-embryonic visceral endoderm, such as APOA2, CITED1 and TTR (Nowotschin et al. 2019).

There was insufficient evidence to annotate the additional 26/87 cell types present in the mouse atlas.

#### Cell type colours

Cell type colours were chosen to most easily facilitate comparisons between similar cell types in the mouse. For 52/67 annotations, where there is a confident or suspected 1-1 mapping between cell types, consistent colours was used. Additional colours were chosen with consideration to visual clarity and palettes used for similar tissues.

### References

See `references.bib`.

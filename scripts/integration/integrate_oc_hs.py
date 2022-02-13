# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 11:46:00 2021

@author: Daniel Keitley
"""

from samap.mapping import SAMAP
from samap.analysis import get_mapping_scores, GenePairFinder
from samalg import SAM
from samap.utils import save_samap, load_samap
import pandas as pd
import scanpy as sc
import numpy as np
import csv


data_path = "../../data-in/integration/"

# Human
h_ensembl = pd.read_csv(data_path + "human_transcript_data.tsv",sep="\t")
h_genes = pd.read_csv(data_path + "human_genes.tsv",sep="\t")
h_names = [(y,x) for x in h_genes["x"] for y in h_ensembl.loc[h_ensembl["external_gene_name"]==x,"ensembl_transcript_id_version"] ]

# Rabbit
with open(data_path +'r_names.tsv') as f:
    r_names=[tuple(line) for line in csv.reader(f,delimiter='\t')]
    

hs_path = data_path + "Tyser_2020/adata.h5ad"
oc_path = "../../data-in/rabbit/anndata.h5ad"

sm = SAMAP(hs_path, oc_path, "hs", "oc",
                   f_maps = data_path + 'maps/',
                              names1=h_names,names2=r_names,
                                         save_processed=False)


sm_out = sm.run()
sm_out.adata.write("../../data-out/integration/oc_hs_samap.h5ad")
save_samap(sm_out,"../../data-out/integration/oc_hs_samap")

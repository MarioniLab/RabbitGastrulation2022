# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 14:41:55 2021

@author: Daniel Keitley
"""

import scanpy as sc
import pandas as pd

def loadRabbitData(data_path="../data-in/rabbit/anndata.h5ad"):
    r_data = sc.read_h5ad(data_path)
    return(r_data)

def loadMouseData(data_path="../data-in/mouse/anndata.h5ad"):
    m_data = sc.read_h5ad(data_path)
    return(m_data)


def exportGMT(adata, group, export):
    for ob in adata.obs[group].cat.categories:
        df = pd.DataFrame([ob,"-"] + adata.obs[adata.obs[group]==ob].index.values.tolist())
        df.T.to_csv(export, mode='a', sep="\t", header=False, index=False)
    
    
def updateAnndata():
    """
    Update Anndata object with new metadata.

    """
    
    

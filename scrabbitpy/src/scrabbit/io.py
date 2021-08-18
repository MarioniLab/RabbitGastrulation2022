# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 14:41:55 2021

@author: Daniel Keitley
"""

import scanpy as sc
import pandas as pd
import anndata

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
    

def loadWOTCSV(wot_path):
    adata = anndata.read_csv(wot_path + "X.csv",delimiter="\t")
    adata.obs = pd.read_csv(wot_path + "obs.csv", sep="\t",index_col=0)
    adata.var = pd.read_csv(wot_path + "var.csv", sep="\t",index_col=0)
    return(adata)
    

def updateAnndata():
    """
    Update Anndata object with new metadata.

    """
    
    

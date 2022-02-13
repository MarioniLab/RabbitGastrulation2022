# -*- coding: utf-8 -*-
"""
Created on Sat Aug 21 17:08:27 2021

@author: Daniel Keitley
"""

import pandas as pd
import numpy as np
import scipy
import scanpy as sc


def calcGeneSpecificity(adata,group_by=None, mask=None):
    if(mask is None):
        mask = pd.get_dummies(adata.obs[group_by])
        
    ncells=mask.sum(0)
    ncells_array = np.squeeze(np.asarray(1/ncells))
    ncells_mat = scipy.sparse.diags(ncells_array)
    mask = scipy.sparse.csr_matrix(mask).T

    ctype_sum = mask@adata.X
    ctype_mean = ncells_mat@ctype_sum
    
    N = ctype_mean.shape[0]
    row_sums = np.squeeze(np.asarray(ctype_mean.sum(axis=1)))[:,None]
    gspec = ctype_mean/(row_sums/N) #TODO: make this sparse    
    
    return(gspec)



# Extract one-to-one orthologs in both filtered datasets
def chooseCommonGenes(r_filt,m_filt,orthologs,use_hvgs=False,join_type="intersect"):
    r_genes = r_filt.var.index
    m_genes = m_filt.var.index
    
    if(use_hvgs):
        r_genes = r_genes[r_filt.var["highly_variable"]]
        m_genes = m_genes[m_filt.var["highly_variable"]]
        
    if(join_type=="union"): 
        genes_in_both = orthologs.loc[orthologs["rabbit"].isin(r_genes) | orthologs["mouse"].isin(m_genes),:]
    else:
        # TODO: Make sure features are in common after this
        genes_in_both = orthologs.loc[orthologs["rabbit"].isin(r_genes) & orthologs["mouse"].isin(m_genes),:]
    
    return(genes_in_both)




def runGeneSpecificity(rabbit,mouse,orthologs,ctype_shared,hvgs=None):
    r_filt = filterGenes(rabbit,orthologs["rabbit"],ctype_shared,hvgs=hvgs)
    m_filt = filterGenes(mouse,orthologs["mouse"],ctype_shared,hvgs=hvgs)
    
    use_hvgs=False
    if hvgs is not None:
        use_hvgs = True
        
    genes_in_both = getCommonOrthologs(r_filt,m_filt,orthologs,use_hvgs=use_hvgs,join_type="intersect")
    r_filt = r_filt[:,genes_in_both["rabbit"]]
    m_filt = m_filt[:,genes_in_both["mouse"]]
    
    r_gspec = calcGeneSpecificity(r_filt,"celltype")
    r_gspec = np.array(r_gspec)
    
    m_gspec = calcGeneSpecificity(m_filt,"celltype")
    m_gspec = np.array(m_gspec)
    
    N = len(ctype_shared)
    gspec_cor = np.corrcoef(r_gspec,m_gspec)
    gspec_cor = gspec_cor[0:N,N:]
    
    ctype_list = np.sort(list(ctype_shared)) # ctype_shared a set so careful of ordering
    gspec_cor = pd.DataFrame(gspec_cor,columns=ctype_list,index=ctype_list)

    return(gspec_cor)
    

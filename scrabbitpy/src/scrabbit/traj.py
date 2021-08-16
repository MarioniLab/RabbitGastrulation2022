# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 16:12:41 2021

@author: Daniel Keitley
"""

import scanpy as sc
import pandas as pd
import wot
import numpy as np

def computeGeneScores(raw_adata, corrected_adata, gs_path,save=False, out_path=""):
    gs = wot.io.read_sets(gs_path, raw_adata.var.index.values)
    gs_scores = pd.DataFrame(index=corrected_adata.obs.index)
    for j in range(gs.shape[1]):
        gene_set_name = str(gs.var.index.values[j])
        result = wot.score_gene_sets(ds=raw_adata, gs=gs[:, [j]], permutations=0, method='mean_z_score')
        gs_scores[gene_set_name] = result['score']
    if(save):
        gs_scores.to_csv(out_path, index_label='id')
    return(gs_scores)


def logistic(x, L, k, x0=0):
    f = L / (1 + np.exp(-k * (x - x0)))
    return f

def gen_logistic(p, beta_max, beta_min, pmax, pmin, center, width):
    return beta_min + logistic(p, L=beta_max - beta_min, k=4 / width, x0=center)

def beta(p, beta_max=1.7, beta_min=0.3, pmax=1.0, pmin=-0.5, center=0.25):
    return gen_logistic(p, beta_max, beta_min, pmax, pmin, center, width=0.5)

def delta(a, delta_max=1.7, delta_min=0.3, amax=0.5, amin=-0.4, center=0.1):
    return gen_logistic(a, delta_max, delta_min, amax, amin, center,
                          width=0.2)

def computeGrowthRates(gene_set_scores,save=False,out_path=""):
    proliferation=gene_set_scores['PROLIFERATION']
    apoptosis = gene_set_scores['CELL_DEATH']
    
    birth = beta(proliferation)
    death = delta(apoptosis)

    # growth rate is given by 
    gr = np.exp(birth-death)
    growth_rates_df = pd.DataFrame(index=gene_set_scores.index, data={'cell_growth_rate':gr})
    
    #if(save):
     #   growth_rates_df.to_csv('data/growth_gs_init.txt')
        
    return(growth_rates_df)



def runWOT(adata, gs_path, out_path, epsilon=0.05, lambda1=1, lambda2=50,
           local_pca=0, growth_iters=3):
    
    # Create AnnData of batch corrected PCs
    corrected = sc.AnnData(adata.obsm["X_pca"],obs=adata.obs)
    
    # Compute gene scores
    g_scores = computeGeneScores(adata, corrected, save=True,
                                 gs_path = gs_path,
                                 out_path = out_path + "growth_gene_set_scores.csv")
    
    # Compute cell growth rates
    g_rates = computeGrowthRates(g_scores, save=True,
                              out_path = out_path + "growth_gs_init.txt")
    corrected.obs["cell_growth_rate"] = g_rates
    
    
    # Create WOT model
    model = wot.ot.OTModel(corrected,epsilon=epsilon,lambda1=lambda1,
                             lambda2=lambda2,local_pca=local_pca,
                             growth_iters=growth_iters)
    
    # Compute transport maps
    model.compute_all_transport_maps(tmap_out = out_path + "tmaps/")


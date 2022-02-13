# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 16:15:20 2021

@author: Daniel Keitley
"""

import scrabbit

# Load data
m_data = scrabbit.io.loadMouseData("../../data-in/mouse/anndata.h5ad")

in_path = "../../data-in/trajectory_analysis/mouse/"
out_path = "../../data-out/trajectory_analysis/mouse/"

scrabbit.traj.runWOT(m_data, in_path + "growth_gene_sets.gmx", out_path)

# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 17:17:21 2021

@author: Daniel Keitley
"""

import scrabbit

# Load data
r_data = scrabbit.io.loadRabbitData("../../data-in/rabbit/anndata.h5ad")

stage_dict = {"GD7":7,"GD8":8,"GD9":9}
r_data.obs['day'] = r_data.obs['stage'].replace(stage_dict)

in_path = "../../data-in/trajectory_analysis/rabbit/"
out_path = "../../data-out/trajectory_analysis/rabbit/"

scrabbit.traj.runWOT(r_data, in_path + "growth_gene_sets.gmx", out_path)

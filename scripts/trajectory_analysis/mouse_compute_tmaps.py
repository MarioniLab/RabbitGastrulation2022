# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 16:15:20 2021

@author: Daniel Keitley
"""

import scrabbit


# Load data
m_data = scrabbit.io.loadMouseData("../../data-in/mouse/anndata.h5ad")

# Index with cell names
m_data.obs.index = m_data.obs["cell"]

# TODO: Add to default anndata
# Add time information
stage_dict = {"E6.5":6.5,"E6.75":6.75,"E7.0":7,"E7.25":7.25,"E7.5":7.5,
                    "E7.75":7.75,"E8.0":8,"E8.25":8.25,"E8.5":8.5,
                    "mixed_gastrulation":None,"E8.75":8.75,"E9.0":9,
                    "E9.25":9.25,"E9.5":9.5}

m_data.obs['day'] = m_data.obs['stage'].replace(stage_dict)
m_data.obs['day'] = m_data.obs['day'].astype('float') # needed for 0.5 days


in_path = "../../data-in/trajectory_analysis/mouse/"
out_path = "../../data-out/trajectory_analysis/mouse/"

scrabbit.traj.runWOT(m_data, in_path + "growth_gene_sets.gmx", out_path)

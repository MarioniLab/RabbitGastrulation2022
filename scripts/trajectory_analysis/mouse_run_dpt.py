# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 16:54:02 2021

@author: Daniel Keitley
"""


import scrabbit

out_path = "data-out/trajectory_analysis/mouse/dpt/"
Path(out_path).mkdir(parents=True, exist_ok=True)

m_data = scrabbit.io.loadMouseData("data-in/mouse/anndata.h5ad")

m_data = scrabbit.traj.runDPT(m_data, "celltype", "Epiblast", denoise=True,
                     n_neighbours=100)

m_data.write(out_path + "anndata.h5ad")
m_data.obs["dpt_pseudotime"].to_csv(out_path + "dpt_pseudotime.tsv",sep="\t")
m_data.obs["dpt_groups"].to_csv(out_path + "dpt_groups.tsv",sep="\t")

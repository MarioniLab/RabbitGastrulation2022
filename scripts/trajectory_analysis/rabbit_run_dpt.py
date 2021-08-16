# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 16:48:42 2021

@author: Daniel Keitley
"""


import scrabbit

out_path = "data-out/trajectory_analysis/rabbit/dpt/"
Path(out_path).mkdir(parents=True, exist_ok=True)

r_data = scrabbit.io.loadRabbitData("data-in/rabbit/anndata.h5ad")

r_data = scrabbit.traj.runDPT(r_data, "celltype", "Epiblast", denoise=True,
                     n_neighbours=100)

r_data.write(out_path + "anndata.h5ad")
r_data.obs["dpt_pseudotime"].to_csv(out_path + "dpt_pseudotime.tsv",sep="\t")
r_data.obs["dpt_groups"].to_csv(out_path + "dpt_groups.tsv",sep="\t")
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 18:18:14 2021

@author: Daniel Keitley
"""
import scrabbit
import os.path
from pathlib import Path

out_path = "data-out/trajectory_analysis/mouse/wot/"
Path(out_path).mkdir(parents=True, exist_ok=True)

celltypes_path = "data-in/trajectory_analysis/mouse/celltypes.gmt"

# Check cell types GMT file exists
if(~os.path.isfile(celltypes_path)):
    m_data = scrabbit.io.loadMouseData("data-in/mouse/anndata.h5ad")
    scrabbit.io.exportGMT(m_data, "celltype", export=celltypes_path)   


scrabbit.traj.computeWOTProbabilities(tmap_path = "data-out/trajectory_analysis/mouse/tmaps/",
                                      celltypes_path = celltypes_path,
                                      export_dir = out_path)

# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 17:33:00 2021

@author: Daniel Keitley
"""

import scrabbit
import os.path

out_path = "data-out/trajectory_analysis/rabbit/wot/"
Path(out_path).mkdir(parents=True, exist_ok=True)

celltypes_path = "data-in/trajectory_analysis/rabbit/celltypes.gmt"

# Check cell types GMT file exists
if(~os.path.isfile(celltypes_path)):
    r_data = scrabbit.io.loadRabbitData("data-in/rabbit/anndata.h5ad")
    scrabbit.io.exportGMT(r_data, "celltype", export=celltypes_path)   


scrabbit.traj.computeWOTProbabilities(tmap_path = "data-out/trajectory_analysis/rabbit/tmaps/",
                                      celltypes_path = celltypes_path,
                                      export_dir = out_path)

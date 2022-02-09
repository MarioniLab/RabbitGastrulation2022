# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 17:33:00 2021

@author: Daniel Keitley
"""

import scrabbit
from pathlib import Path

out_path = "data-out/trajectory_analysis/macaque/wot/"
Path(out_path).mkdir(parents=True, exist_ok=True)

celltypes_path = "data-in/trajectory_analysis/macaque/celltypes.gmt"

scrabbit.traj.computeWOTProbabilities(tmap_path = "data-out/trajectory_analysis/macaque/tmaps/",
                                      celltypes_path = celltypes_path,
                                      end_timepoint = 14,
                                      export_dir = out_path)

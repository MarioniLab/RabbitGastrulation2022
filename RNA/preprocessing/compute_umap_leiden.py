import scanpy as sc
import pandas as pd
import numpy as np

seed = 64

# Load data
r_data = sc.read_h5ad("data-in/r_data_processed.h5ad")

r_pcs = pd.read_csv("data-in/corrected_pcs.tsv", sep="\t")
r_data.obsm["X_pca"] = np.array(r_pcs)

# Compute UMAP
sc.pp.neighbors(r_data, n_neighbors=300, use_rep="X_pca", random_state=seed)
sc.tl.umap(r_data, min_dist=0.9, random_state=seed)

# Cluster
sc.tl.leiden(rabbit,resolution=1, key_added="leiden_res1")
sc.tl.leiden(rabbit,resolution=2, key_added="leiden_res2")
sc.tl.leiden(rabbit,resolution=1.5, key_added="leiden_res1_5")
sc.tl.leiden(rabbit,resolution=2.5, key_added="leiden_res2_5")
sc.tl.leiden(rabbit,resolution=3, key_added="leiden_res3")
sc.tl.leiden(rabbit,resolution=5, key_added="leiden_res5")
sc.tl.leiden(r_data,resolution=6, key_added="leiden_res6")
sc.tl.leiden(r_data,resolution=7, key_added="leiden_res7")
sc.tl.leiden(r_data,resolution=8, key_added="leiden_res8")
sc.tl.leiden(rabbit,resolution=10,key_added="leiden_res10")

# Export
np.savetxt("r_umap.tsv", r_data.obsm["X_umap"], delimiter="\t")
r_data.write("data-out/r_umap_leiden.h5ad")


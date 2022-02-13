import scrabbit
import scanpy as sc

# Load data
mf_data = sc.read_h5ad("../../../data-in/datasets/Yang_2021/anndata.h5ad")

stage_dict = {"d10":10, "d12":12, "d14":14}

mf_data.obs['day'] = mf_data.obs['stage'].replace(stage_dict)

in_path = "../../../data-in/trajectory_analysis/macaque/"
out_path = "../../../data-out/trajectory_analysis/macaque/"

scrabbit.traj.runWOT(mf_data, in_path + "growth_gene_sets.gmx", out_path)

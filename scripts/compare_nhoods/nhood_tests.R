

library(optparse)
library(scrabbitr)

set.seed(44)


###################
## Parse options ##
###################

option_list = list(
  make_option(c("-n1", "--milo1"), type="character", default=NULL,
              help="path to first milo file", metavar="character"),
  make_option(c("-n2", "--milo2"), type="character", default=NULL,
              help="path to second milo file", metavar="character"),
  make_option(c("-t", "--test"), type="character", default=NULL,
              help="name of the test to run", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="All",
              help="basic outpath", metavar="character"))

opt_parser = OptionParser(option_list=option_list)
opts = parse_args(opt_parser)



# Load data
r_milo <- readRDS(opts$milo1)
m_milo <- readRDS(opts$milo2)



runNhoodTest(opts$test, exp_outdir,
                         r_milo, m_milo, orthologs, ...)



# hvg_selection = "scran" | "seurat" | (provide file of genes)
# hvg_exclude = (provide file of genes)
# hvg_join_type = "union" | "intersection"
# scaling = TRUE | FALSE
# sim_preprocessing = "gene_spec" | NULL
# sim_measure = "pearson" | "spearman" | "cosine" | "euclidean" | "mnn"




gut_celltypes <- c("Visceral endoderm","ExE endoderm", "Midgut", "Gut tube","Primitive Streak", 
                   "Epiblast","Foregut", "Hindgut","Definitive Endoderm","Anterior Primitive Streak",
                   "Hypoblast","Visceral YS endoderm","Gut endoderm","Def. Endoderm","Ectoderm")

meso_celltypes <- c("Sclerotome","Dermomyotome", "Somitic mesoderm", "Anterior somitic tissues",
                    "Posterior somitic tissues","Cranial mesoderm","Anterior Primitive Streak",
                    "Caudal epiblast","Paraxial mesoderm","Intermediate mesoderm","Epiblast",
                    "Primitive Streak","Presomitic mesoderm","NMPs","NMPs/Mesoderm-biased",
                    "Caudal mesoderm", "Nascent mesoderm")

allan_celltypes <- c("Allantois","Lateral plate mesoderm", "Epiblast", "Primitive Streak", "Nascent mesoderm")




exportTrajectoryTestPlots <- function(r_milo, m_milo, df_simFilt, export_dir) {
  
  # Gut
  plot_gut <- plotTrajMappings(r_milo, m_milo, df_simFilt, 
                               group="celltype", groups=gut_celltypes, 
                               dimred="UMAP",
                               offset=c(0,3),line_alpha=0.15,edge_alpha=0.01, 
                               reflect.X=FALSE, legend_pos="right")
  
  ggsave(paste0(export_dir,"r_gut_sim.pdf"), plot_gut, width=18, 
         height=8, dpi=300)
  
  
  # Mesoderm
  plot_meso <- plotTrajMappings(r_milo, m_milo, df_simFilt, 
                                group="celltype", groups=meso_celltypes, 
                                dimred="UMAP", 
                                offset=c(20,20),reflect.X=TRUE,line_alpha=0.1,
                                edge_alpha=0.01,legend_pos="right")
  
  ggsave(paste0(exp_outdir, "r_mesoderm_sim.pdf"), plot_meso, width=18, 
         height=8, dpi=300)
  
  
  # Allantois
  plot_allan <- plotTrajMappings(r_milo, m_milo, df_simFilt, 
                                 group="celltype", groups=allan_celltypes, 
                                 dimred="UMAP", offset=c(0.5,0),
                                 reflect.X=F,reflect.Y=T,line_alpha=0.05,
                                 edge_alpha=0.01, legend="right")
  
  ggsave(paste0(exp_outdir, "r_allantois_sim.pdf"), plot_allan, width=18, 
         height=8, dpi=300)
  
  
}


exporNhoodSimTestPlots <- function(r_milo, m_milo, nhood_sim, export_dir) {
  
  
  r_maxNhoods <- getMaxMappings(nhood_sim, 1, long_format=FALSE)
  m_maxNhoods <- getMaxMappings(nhood_sim, 2, long_format=FALSE)
  df_simFilt <- rbind(r_maxNhoods, m_maxNhoods)
  
  
  # Maximum similarity plots
  plot_r_maxsim <- plotNhoodMaxSim(r_milo, r_maxNhoods)
  ggsave(paste0(exp_outdir,"r_max_sim_plot.pdf"), plot_r_maxsim, 
         width=10, height=8, dpi=300)
  
  plot_m_maxsim <- plotNhoodMaxSim(m_milo, m_maxNhoods)
  ggsave(paste0(exp_outdir,"m_max_sim_plot.pdf"), plot_m_maxsim, 
         width=10, height=8, dpi=300)
  
  plot_rm_maxsim <- grid.arrange(p1,p2,nrow=1)
  ggsave(paste0(exp_outdir,"rm_max_sim_plot.pdf"), plot_rm_maxsim, 
         width=18, height=8, dpi=300)
  
  
  # Rideline plot
  plot_ridgeline <- plotNhoodSimGroups(r_milo, r_maxNhoods$sim, xlabel="Correlation", 
                                       ylabel="Cell type",
                                       size=0.1,rel_min_height=0.001)
  ggsave(paste0(exp_outdir,"r_celltype_ridgeline_plot.pdf"), plot_ridgeline, width=6, 
         height=6, dpi=300)
  
  
  # Trajectory plots
  exportTrajectoryTestPlots(r_milo, m_milo, df_simFilt, export_dir)
  
  
}




runNhoodPipeline <- function(export_dir, r_milo, m_milo, rm_orthologs, ...) {
  
  # Create output directory if necessary
  dir.create(file.path(export_dir), showWarnings = FALSE)
  
  out <- calcNhoodSim(r_milo, m_milo, rm_orthologs,
                      export_dir = export_dir, ...)
  
  exporNhoodSimTestPlots(r_milo, m_milo, out$nhood_sim, export_dir)
  
}



# Test effect of number of genes
runNGenesTest <- function(test_dir, r_milo, m_milo, rm_orthologs) {
  
  # Create test directory
  dir.create(file.path(test_dir), showWarnings = FALSE)
  
  n_hvgs <- c(50,100,200,500, 1000, 1500, 2000, 3000, 5000)
  
  for(n_hvg in n_hvgs) {
    test_run_path <- paste0(test_dir, "/", n_hvgs, "_hvgs/")
    runNhoodPipeline(test_run_path, r_milo, m_milo, rm_orthologs, max_hvgs=n_hvg)
  }
  
}




runGSpecTest <- function(test_dir, r_milo, m_milo, rm_orthologs) {
  
  # Create test directory
  dir.create(file.path(test_dir), showWarnings = FALSE)
  
  # Test with/without gene specificity
  runNhoodPipeline(paste0(test_dir, "wout_gspec"), r_milo, m_milo, rm_orthologs, 
                   sim_preprocessing=NULL)
  runNhoodPipeline(paste0(test_dir, "w_gspec"), r_milo, m_milo, rm_orthologs, 
                   sim_preprocessing="gene_spec")
  
}




runNhoodTest <- function(test_name, base_dir,
                            r_milo, m_milo, rm_orthologs, ...) {

  # Create test directory if necessary
  dir.create(file.path(base_dir, test_name), showWarnings = FALSE)
  test_outdir <- paste0(base_dir, test_name, "/")
  
  switch(test_name,
         
         TEST_NGENES = {
           runNGenesTest(test_outdir, r_milo, m_milo, rm_orthologs)
         },
         
         TEST_GSPEC = {
           runGSpecTest(test_outdir, r_milo, m_milo, rm_orthologs)
         },
         
         TEST_CORTYPE = {
           
         },
         
         TEST_HVG_JOIN = {
           
         }
         
  ) 
  

}












# Test with/without gene specificity
runNhoodSimTest("exp-gene_spec", "test-wout_gspec", r_milo, m_milo, rm_orthologs, sim_preprocessing=NULL)
runNhoodSimTest("exp-gene_spec", "test-w_gspec", r_milo, m_milo, rm_orthologs, sim_preprocessing="gene_spec")


# Test pearson vs spearman
runNhoodSimTest("exp-cor_type", "test_pearson", r_milo, m_milo, rm_orthologs, sim_measure="pearson")
runNhoodSimTest("exp-cor_type", "test_spearman", r_milo, m_milo, rm_orthologs, sim_measure="spearman")




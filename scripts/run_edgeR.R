
library(optparse)
library(scrabbitr)

set.seed(44)


###################
## Parse options ##
###################

option_list = list(
  make_option(c("-s", "--sce"), type="character", default=NULL,
              help="path to the sce file", metavar="character"),
  make_option(c("-g", "--groupby"), type="character", default=NULL,
              help="name of grouping variable", metavar="character"),
  make_option(c("-a", "--group1"), type="character", default=NULL,
              help="name of cluster 1", metavar="character"),
  make_option(c("-b", "--group2"), type="character", default="All",
              help="name of cluster 2 (default: All)", metavar="character"),
  make_option(c("-k", "--block"), type="character", default="All",
              help="name of blocking factor", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="All",
              help="basic outpath", metavar="character"))

opt_parser = OptionParser(option_list=option_list)
opts = parse_args(opt_parser)


# Load data
sce <- readRDS(opts$sce)

sce.edgeR <- scrabbitr::prepareEdgeR(sce, opts$groupby, 
                                     opts$group1, opts$group2, 
                                     opts$block)

out <- scrabbitr::runEdgeR(sce.edgeR,
                min_detection_rate_per_group=0.1,
                min.logFC=2,
                threshold_fdr=0.1)


# Save results
saveRDS(out, opts$out)



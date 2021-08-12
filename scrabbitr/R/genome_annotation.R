

#' Obtains reads within a window of a given GRange
#'
#' Constrains the reads considered to within a specified window of a GRange.
#' @param bam_file Rsamtools BamFile object
#' @param gene_range GRange object of gene of interest
#' @param offset Offset from the start and end of the specified gene in which to
#'  consider reads
#' @export
loadReadsAroundGRange <- function(bam_file, gene_range, offset=1e+06) {
  read_range <- gene_range+offset
  scan_params <- Rsamtools::ScanBamParam(which = read_range, what = scanBamWhat())
  reads <- GenomicAlignments::readGappedReads(bam_file,param=scan_params)
  return(as(reads, "GRanges"))
}


#' Plots histogram of distances of reads that are nearest to the gene given
#' @importFrom ggplot2 ggplot geom_histogram xlab ggtitle aes
#' scale_x_continuous scale_y_continuous ylab
#' @importFrom IRanges distanceToNearest
#' @export
plotReadDistsNearGene <- function(gene_id, read_ranges, gene_ranges,
                                  colour="seagreen4") {
  read_dists <- distanceToNearest(read_ranges@ranges,gene_ranges@ranges)
  gene_index <- which(gene_ranges$gene_id==gene_id)
  gene_dists <- read_dists[subjectHits(read_dists)==gene_index]
  gene_dists <- as.data.frame(gene_dists)

  p <- ggplot(gene_dists, aes(x=distance/(1*10^3))) +
    geom_histogram(colour=colour,fill=colour,alpha=0.8) + theme_light() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab("Distance of intergenic reads to gene (kb)") +
    ylab("Count") +
    ggtitle(gene_id)
  return(p)
}

# plotReadDistsNearGene <- function(gene_id, bam_file, gtf_gf) {
#   genes <- loadGenesNearGene(gtf_gf, gene_id)
#   reads <- loadReadsAroundGRange(bam_file, gene_id,genes)
#   p = plotReadDistsNearGene(gene_id, reads,genes)
#   return(p)
# }




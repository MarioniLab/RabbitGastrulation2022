
library(here)
library(Rsamtools)
library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(ggplot2)
library(Rgb)
library(biomaRt)
library(data.table)


here::set_here(getwd())

# File paths - modify where appropriate
# Uses 'here' R package to simplify relative paths
# Assumes working directory is source file location 

sample_bam_path <- here("data-in","SIGAC11_possorted_genome_bam.bam")
sample_bai_path <- here("data-in","SIGAC11_possorted_genome_bam.bam.bai")

intergenic_bam_path <- here("data-in","SIGAC11_intergenic_reads.bam")
intergenic_bai_path <- here("data-in","SIGAC11_intergenic_reads.bam.bai")

gtf_path <- here("data-in","Oryctolagus_cuniculus.OryCun2.0.99.gtf")
alignments_path <- here("data-in","hs_alignment_results.tsv") 

out_gtf_ext <- here("data-out","extended_three_prime_annotation.gtf")
out_gtf_align <- here("data-out","extension_plus_alignments.gtf")


# Load BAM files
sample_bam <- BamFile(file=sample_bam_path,
                   index=sample_bai_path)

intergenic_bam <- BamFile(file=intergenic_bam_path,
                   index=intergenic_bai_path)

chr_lengths <- seqlengths(seqinfo(sample_bam))



# Load GTF file

# Loads using both GenomicFeatures package and Rgb package. 
# GenomicFeatures (gf) conveniently loads gene coordinates as GRange objects useful for 
# calculaing gene distances 

# Rgb (rgb) loads all information from GTF file directly including gene_name which seems to
# be missed by GenomicFeatures 

# TODO: Modify functions below to rely on one gtf format 

gtf_gf <- makeTxDbFromGFF(gtf_path, format="gtf")
gtf_genes_gf <- genes(gtf_gf)


gtf_rgb <- read.gtf(gtf_path)
gtf_genes_rgb <- gtf_rgb[(gtf_rgb$feature=="gene"),
                                c("gene_id","gene_name","seqname","strand","start","end")]
gtf_genes_rgb <- as.data.frame(gtf_genes_rgb)
gene_ranges_rgb <- gtf_genes_gf[as.character(gtf_genes_rgb$gene_id),]





# Exploring read distances --------------------------------------------------

# Example queries:

# getReadDistsNearGene("ENSOCUG00000007368",intergenic_bam,gtf_gf) + xlim(0,1) + ggtitle("FOXA2")
# getReadDistsNearGene("ENSOCUG00000017248",intergenic_bam,gtf_gf) + ggtitle("TAL1") + xlim(0,1)
# getReadDistsNearGene("ENSOCUG00000011639",intergenic_bam,gtf_gf) + ggtitle("SOX9") + xlim(0,1)
# getReadDistsNearGene("ENSOCUG00000025597",intergenic_bam,gtf_gf) + ggtitle("GATA1") + xlim(0,1) + ylim(0,25)
# getReadDistsNearGene("ENSOCUG00000023670",intergenic_bam,gtf_gf) + ggtitle("RAMP2") + xlim(0,1) + ylim(0,100) 


# grange <- getGeneRange("ENSOCUG00000007368",gtf_genes_gf) + 1e+04
# getReadDistsWithinRange(grange,intergenic_bam,gtf_gf) + xlim(0,1) + ylim(0,20000)

# grange <- GRanges(seqnames = "4",ranges = IRanges(0,1e+07))
# getReadDistsWithinRange(grange,intergenic_bam,gtf_gf) + xlim(0,1) + ylim(0,20000)


# Used for publication figure
# grange <- GRanges(seqnames = "1",ranges = IRanges(0,1e+07))
# p <- getReadDistsWithinRange(grange,intergenic_bam,gtf_gf)
#
# p + xlim(0,1) + ylim(0,5000) + 
#   geom_histogram(color = "#000000",fill = "royalblue2", alpha=0.6)  + 
#   geom_vline(aes(xintercept = 0.6), color = "#000000", linetype = "dashed", 
#              size = 1)  + expand_limits(x = 0, y = 0) + 
#   scale_y_continuous(expand = c(0, 0), limits=c(0,5000)) + 
#   scale_x_continuous(expand = c(0, 0), limits=c(0,1)) + 
#   theme_bw() + ylab("count") + theme(aspect.ratio = 1)


# Loads gene ranges of all genes on same chromosome as a given gene
loadGenesNearGene <- function(gtf_gf,gene_id) {
  gtf_genes <- genes(gtf_gf)
  gene_index <- which(gtf_genes$gene_id==gene_id)
  gene_chr <- runValue(seqnames(gtf_genes)[gene_index])
  chr_genes = gtf_genes[seqnames(gtf_genes) == gene_chr]
  return(chr_genes)
}


# Gets GRange of a given gene
getGeneRange <- function(gene_id,gene_ranges) {
  grange_data <- gene_ranges[gene_ranges$gene_id == gene_id,]
  gene_range <- ranges(grange_data)
  return(GRanges(seqnames = seqnames(grange_data),
                 ranges=gene_range))
}


# Loads reads within a 1Mb window of a given gene - reduces runtime
loadReadsNearGene <- function(bam_file,gene_id,gene_ranges,offset=1e+06) {
  grange_data <- gene_ranges[gene_ranges$gene_id == gene_id,]
  gene_chr <- seqnames(grange_data)
  gene_range <- ranges(grange_data)
  read_range <- gene_range+offset
  
  scan_range <- GRanges(seqnames = gene_chr,
                        ranges = read_range)
  scan_params <- ScanBamParam(which = scan_range, what = scanBamWhat())
  
  reads <- readGappedReads(bam_file,param=scan_params)
  return(as(reads, "GRanges"))
}


# Plots histogram of distances of reads whose nearest to the gene given
plotReadDistsNearGene <- function(gene_id,read_ranges,gene_ranges) {
  read_dists <- distanceToNearest(read_ranges@ranges,gene_ranges@ranges)
  
  gene_index <- which(gene_ranges$gene_id==gene_id)
  gene_dists <- read_dists[subjectHits(read_dists)==gene_index]
  gene_dists <- as.data.frame(gene_dists)
  p = ggplot(gene_dists, aes(x=distance/(1*10^3))) + 
    geom_histogram() + xlab("Distance of intergenic reads to gene (kb)") + ggtitle(gene_id)
  return(p)
}


getReadDistsNearGene <- function(gene_id,bam_file,gtf_gf) {
  genes <- loadGenesNearGene(gtf_gf,gene_id)
  reads <- loadReadsNearGene(bam_file,gene_id,genes)
  p = plotReadDistsNearGene(gene_id,reads,genes)
  return(p)
}


#TODO: Rename this function - misleading
# Loads all genes on the same chromosome as the grange
loadGenesWithinRange <- function(grange,gtf_gf) {
  genes <- genes(gtf_gf)
  genes <- genes[seqnames(genes)==as.character(runValue(seqnames(grange))),]
  return(genes)
}


loadReadsWithinRange <- function(grange,bam_file) {
  scan_params <- ScanBamParam(which = grange, what = scanBamWhat())
  
  reads <- readGappedReads(bam_file,param=scan_params)
  return(as(reads, "GRanges"))
}


plotReadDistsWithinRange <- function(read_ranges,gene_ranges) {
  read_dists <- distanceToNearest(read_ranges@ranges,gene_ranges@ranges)
  read_dists <- as.data.frame(read_dists)
  p = ggplot(read_dists, aes(x=distance/(1*10^3))) + 
    geom_histogram() + xlab("Distance of intergenic reads to nearest gene (kb)") 
  return(p)
}


getReadDistsWithinRange <- function(grange,bam_file,gtf_gf) {
  genes <- loadGenesWithinRange(grange,gtf_gf)
  reads <- loadReadsWithinRange(grange,bam_file)
  p <- plotReadDistsWithinRange(reads,genes) + ggtitle(as.character(grange))
  return(p)
}




# Adding 3' extension -----------------------------------------------------


createArtificialExon <- function(gene_id,gene_name,gene_seqname,exon_start,exon_end,gene_strand) {
  if(is.na(gene_name)) {
    gene_name = gene_id
  }
  
  new_exon <- c(gene_seqname,"ensembl","exon", exon_start, exon_end,".",gene_strand,".",
                paste0('gene_id "',gene_id,
                       '"; gene_version "3"; transcript_id "ARTIFICIAL_TRANSCRIPT_',gene_id,'"; transcript_version "3"; exon_number "999"; gene_name "',
                       gene_name,'"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "ARTIFICIAL_TRANSCRIPT_',
                       gene_id,'"; transcript_source "ensembl"; transcript_biotype "protein_coding"; exon_id "ARTIFICIAL_EXON_',
                       gene_id,'"; exon_version "2";'))
  return(new_exon)
}



createArtificialCDS <- function(gene_id,gene_name,gene_seqname,exon_start,exon_end,gene_strand) {
  if(is.na(gene_name)) {
    gene_name = gene_id
  }
  new_cds <- c(gene_seqname,"ensembl","CDS", exon_start, exon_end,".",gene_strand,1,
               paste0('gene_id "',
                      gene_id,'"; gene_version "3"; transcript_id "ARTIFICIAL_TRANSCRIPT_',
                      gene_id,'"; transcript_version "3"; exon_number "999"; gene_name "',
                      gene_name,'"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "ARTIFICIAL_TRANSCRIPT_',
                      gene_id,'"; transcript_source "ensembl"; transcript_biotype "protein_coding"; protein_id "ARTIFICIAL_PROTEIN_',
                      gene_id,'"; protein_version "1";'))
  return(new_cds)
  
}


writeToGtf <- function(gtf_path,entry) {
  write.table(entry, file = gtf_path, 
              append = TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}



# Adds artifical exons to the 3' end of genes 
# Expects 'genes' in gtf_genes_rgb format i.e. gene_id, gene_name,...
# gene_ranges in format given by gtf_genes_gf
# extend_dist must be given in kbp
# e.g. addArtificialExons(gtf_genes_rgb, gene_ranges_rgb, chr_lengths, out_gtf_ext)
addArtificialExons <- function(genes, gene_ranges, chr_lengths, new_gtf_path, extend_dist=0.6) {

  apply(genes,1,function(gene) {
    gene_id <- gene[["gene_id"]]
    gene_name <- gene[["gene_name"]]
    gene_strand <- gene[["strand"]]
    gene_seqname <- gene[["seqname"]]
    gene_range <- gene_ranges[as.factor(gene_id),]
    
    start_of_gene <- as.numeric(as.data.frame(gene_range@ranges)[["start"]])
    end_of_gene <- as.numeric(as.data.frame(gene_range@ranges)[["end"]])
    
    
    if(as.character(gene_strand)=="+") {
      new_exon_start <- end_of_gene
      new_exon_end <- min(end_of_gene + extend_dist*1000,chr_lengths[[gene_seqname]])
      next_gene <- precede(gene_range,gene_ranges,ignore.strand=FALSE)
    } else {
      new_exon_start <- max(start_of_gene - extend_dist*1000,1)
      new_exon_end <- start_of_gene
      next_gene <- follow(gene_range,gene_ranges,ignore.strand=FALSE) 
    }
    
    
    # Last gene on chromosome
    if(is.na(next_gene)) {
      dist_next_gene <- chr_lengths[[gene_seqname]] - end_of_gene
    } else {
      # Note ignore.strand=FALSE means added exons can technically overlap. 
      dist_next_gene <- distance(gene_range,gene_ranges[next_gene,],ignore.strand=FALSE)
    }
    
    # Check intergenic distance is greater than twice extension the length
    if((dist_next_gene > 2*extend_dist*1000)) {
      
      new_exon <- createArtificialExon(gene_id,gene_name,gene_seqname,new_exon_start,new_exon_end,gene_strand)
      new_cds <- createArtificialCDS(gene_id,gene_name,gene_seqname,new_exon_start,new_exon_end,gene_strand)
      
      writeToGtf(new_gtf_path,t(new_exon))
      writeToGtf(new_gtf_path,t(new_cds))

    }

  })
  
}



# Adding aligned exons -------------------------------------------------------


r_alignments <- read.table(alignments_path,stringsAsFactors = F,header = TRUE)

createArtificialGene <- function(gene_id,gene_name,gene_seqname,gene_start,gene_end,gene_strand) {
  if(is.na(gene_name)) {
    gene_name = gene_id
  }
  
  new_gene <- c(gene_seqname,"ensembl","gene", gene_start, gene_end,".",gene_strand,".",
                paste0('gene_id "',gene_id, '"; gene_version "3"; gene_name "',
                       gene_name,'"; gene_source "ensembl"; gene_biotype "protein_coding";'))
  return(new_gene)
  
}


mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


getNonOverlappingIntervals <- function(seq_starts,seq_ends) {
  
  sorted_seqs <- sort(seq_starts,index.return=TRUE)
  seq_starts <- sorted_seqs$x
  seq_ends <- seq_ends[sorted_seqs$ix]
  
  cur_exon <- 1
  cur_exon_start <- seq_starts[cur_exon]
  cur_exon_end <- seq_ends[cur_exon]

  new_exon_starts <- c()
  new_exon_ends <- c()
  
  while(cur_exon_end!=-1) {
    next_exon <- cur_exon+1
    next_exon_start <- seq_starts[next_exon]
    next_exon_end <- seq_ends[next_exon]
    
    if(cur_exon_end == max(seq_ends)) {
      new_exon_starts = c(new_exon_starts,cur_exon_start)
      new_exon_ends = c(new_exon_ends,cur_exon_end)
      cur_exon_end=-1
      break()
    }
    
    if(next_exon_start <= cur_exon_end) {
      cur_exon_end <- max(cur_exon_end,next_exon_end)
    } else {
      new_exon_starts <- c(new_exon_starts,cur_exon_start)
      new_exon_ends <- c(new_exon_ends,cur_exon_end)
      cur_exon_start <- next_exon_start
      cur_exon_end <- next_exon_end
    }
    
    cur_exon <- next_exon
    
  }
  
  
  return(data.frame(new_exon_starts,new_exon_ends))
  
}


ext_gtf_granges <- makeGRangesFromDataFrame(ext_gtf,keep.extra.columns = TRUE)
ext_gtf_exons <- gtf_granges[gtf_granges$feature=="exon",]
ext_gtf_genes <- gtf_granges[gtf_granges$feature=="gene",]
out<-addUnannotatedExons(out_gtf_align, r_alignments, ext_gtf_genes, ext_gtf_exons)

# Adds artificial exon to coordinates given by alignment coordinates in new_exons
# Example - out<-addUnannotatedExons(out_gtf_align,rabbit_alignments,gtf_genes_gf)
# Expects columns of r_alignments to correspond with example in alignments.tsv
addUnannotatedExons <- function(new_gtf_path,new_exons,gene_ranges,exon_ranges) {
  exons_added <- 0
  sapply(unique(new_exons[[1]]),function(gene_id) {
    
    new_exons <- unique(new_exons[new_exons[[1]]==gene_id,])
    
    gene_strand <- mode(new_exons$rabbit_align_strand)
    if(gene_strand==1) {gene_strand<-"+"} else {gene_strand<-"-"}
    
    gene_start <- min(new_exons$rabbit_align_start)
    gene_end <- max(new_exons$rabbit_align_end)
    gene_name <- new_exons$human_gene_name[1]
    
    
    # check all exons on same chromsome
    if(length(unique(new_exons$rabbit_align_seq))!=1) {
      stop(paste0("Exons for the same gene mapping to different chromosomes. (Gene ID: ",gene_id,")"))
    }
    
    # check new gene doesn't overlap existing annotation
    gene_seqname <- new_exons$rabbit_align_seq[1]
    grange <- GRanges(seqname=gene_seqname,IRanges(start=gene_start,end=gene_end),strand=gene_strand)
    gene_hits <-distanceToNearest(grange,gene_ranges)
    if(length(gene_hits)==0) {dist_to_gene<-Inf}
    else {dist_to_gene <- elementMetadata(gene_hits)$distance}
  
    if(dist_to_gene!=0) {
      new_gene <- createArtificialGene(gene_id,gene_name,gene_seqname, gene_start,gene_end,gene_strand)
      writeToGtf(new_gtf_path,t(new_gene))
    } else {
      print(paste0("New gene annotation (", gene_id,") overlaps exisiting annotation. Will try annotating exons only..."))
    }
    
    final_exons <- getNonOverlappingIntervals(new_exons$rabbit_align_start,new_exons$rabbit_align_end)
    
    
    for(i in 1:nrow(final_exons)) {
      exon <- final_exons[i,]
      
      exon_start <- exon[[1]]
      exon_end <- exon[[2]]
      
      grange <- GRanges(seqname=gene_seqname,IRanges(start=exon_start,end=exon_end),strand=gene_strand)
      exon_hits <-distanceToNearest(grange,exon_ranges)
      if(length(exon_hits)==0) {dist_to_exon<-Inf}
      else {dist_to_exon <- elementMetadata(exon_hits)$distance}
      
      if(dist_to_exon==0) {
        print("New exon overlaps existing annotation. Skipping exon...")
        print(paste0("Gene ID: ", gene_id,", Seq: ", gene_seqname, ", Start: ",exon_start,", End: ",exon_end))
        next
      }

      new_exon <- createArtificialExon(gene_id,gene_name,gene_seqname,exon_start,exon_end,gene_strand)
      new_cds <- createArtificialCDS(gene_id,gene_name,gene_seqname,exon_start,exon_end,gene_strand)
      
      writeToGtf(new_gtf_path,t(new_exon))
      writeToGtf(new_gtf_path,t(new_cds))
      exons_added <- exons_added + 1
      
    }
    
  })
  
  print(paste0(exons_added, " exons added to annotation. "))
}




# Calculate mapping improvements ------------------------------------------


gr <- GRanges(seqnames = "2",
              ranges = IRanges(start = 1, end = chr_lengths[2]))
params <- ScanBamParam(which = gr, what = scanBamWhat())
reads_sample <- scanBam(sample_bam,param=params)
reads_sample[[1]]$end = reads_sample[[1]]$pos + reads_sample[[1]]$qwidth
read_ranges <- makeGRangesFromDataFrame(reads_sample[[1]],start.field = "pos",end.field = "end",seqnames.field = "rname" )

orig_gtf <- gtf_rgb
ext_gtf <- read.gtf(here("data-out","extended_three_prime_annotation.gtf"))
ext_aln_gtf <- read.gtf(here("data-out","extension_plus_alignments.gtf"))


# e.g. calcPercReadsMappedToGTF(orig_gtf,read_ranges,gr)
calcPercReadsMappedToGTF <- function(gtf,read_ranges,gr) {
  gtf_granges <- makeGRangesFromDataFrame(gtf,keep.extra.columns = TRUE)
  gtf_exons <- gtf_granges[gtf_granges$feature=="exon",]
  gtf_exons <- subsetByOverlaps(gtf_exons,gr)
  
  exon_dists <- distanceToNearest(read_ranges,gtf_exons)
  exon_reads <- exon_dists[elementMetadata(exon_dists)[,1]==0,]
  
  total_reads <- nrow(as.data.frame(read_ranges))
  perc_mapped <- length(queryHits(exon_reads))/total_reads
  
  return(perc_mapped)
  
}






# Exploring read coverage -------------------------------------------------

read_coverage_path <- here("data-in","SIGAC11_coverage.tsv")

con = file(read_coverage_path, "r")
chunk_size <- 100000

read = TRUE
while(read) {
  read_coverage <- read.table(con,header=FALSE,sep="\t",nrows = chunk_size)
}

read_coverage <- fread(read_coverage_path,header = FALSE,sep = "\t",nrows = chunk_size)
colnames(read_coverage) <- c("seqname","position","coverage")

read_coverage_grange <- makeGRangesFromDataFrame(read_coverage,keep.extra.columns = TRUE,
                                            start.field = "position",end.field="position",
                                            seqnames.field="seqname")

complete_coverage <- rep(0,max(read_coverage$position))
complete_coverage[read_coverage$position] = read_coverage$coverage

ext_aln_gtf <- read.gtf(here("../out","extension_plus_alignments.gtf"))
gtf_granges <- makeGRangesFromDataFrame(ext_aln_gtf,keep.extra.columns = TRUE)


non_gtf_coverage <- read_coverage_grange[!read_coverage_grange %over% gtf_granges,]
non_gtf_coverage <- as.data.frame(intg_coverage)
#ggplot(intg_coverage_df,aes(x=start,y=coverage)) + geom_point()
#addReadClusterAnnotations(complete_coverage,read_coverage_grange,gtf_granges)

addReadClusterAnnotations <- function(complete_coverage,read_coverage,gtf_granges,
                                      min_coverage=250,gap=100,cutoff=0.1) {
  
  non_gtf_coverage <- read_coverage[!read_coverage %over% gtf_granges,]
  lrcs <- as.data.frame(non_gtf_coverage[non_gtf_coverage$coverage>min_coverage,])
  gaps <- diff(lrcs$start)
  cluster_gaps <- which(gaps>gap)
  
  start_index <- 1
  for(i in 1:(length(cluster_gaps)+1)) {
    
    end_index<-0
    if(i > length(cluster_gaps)) {end_index<-nrow(lrcs)
    } else { end_index <- cluster_gaps[i]}
    

    lrcs_reads<-lrcs[start_index:end_index,]
    
    clust_max <- lrcs_reads[which.max(lrcs_reads$coverage),]
    clust_thresh <- cutoff*clust_max$coverage
    below_thresh <- which(complete_coverage < clust_thresh)
    
    clust_start <- clust_max$start - min(clust_max$start-below_thresh[below_thresh<clust_max$start])
    clust_end <- clust_max$start + min(below_thresh[below_thresh>clust_max$start]-clust_max$start)
    
    print(clust_start)
    print(clust_end)
    # check this doesn't overlap existing annotation
    
    if(i > length(cluster_gaps)) {start_end<- -1
    } else {start_index <- cluster_gaps[i] + 1}
    
  }
  
}



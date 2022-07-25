here::i_am("01_create_arrow.R")
source(here::here("settings.R"))

# I/O
io$output.directory <- file.path(io$basedir,"ArchR")
setwd(io$output.directory)

# Define arguments
# p <- ArgumentParser(description='')
# p$add_argument('--arrow_files',     type="character",  nargs='+',      help='Arrow files')
# args <- p$parse_args(commandArgs(TRUE))

addArchRThreads(threads = 1) 

# Create project
genomeAnnotation = readRDS(file.path(io$basedir, 'genomeAnnotation.rds'))
geneAnnotation = readRDS(file.path(io$basedir, 'geneAnnotation.rds'))

ArrowFiles = file.path(io$output.directory, list.files(io$output.directory, pattern ='arrow'))

proj <- ArchRProject(
  #ArrowFiles = args$arrow_files, 
  ArrowFiles = ArrowFiles, 
  outputDirectory = "Project",
  copyArrows = TRUE, #This is recommened so that you maintain an unaltered copy for later usage.
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation
)

proj <- filterDoublets(ArchRProj = proj)

saveArchRProject(ArchRProj = proj)

write.table('', file=file.path(io$output.directory,'02_completed.txt'))
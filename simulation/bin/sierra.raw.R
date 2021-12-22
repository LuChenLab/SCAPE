# setwd(
#   "/mnt/data8/zhouran/proj/2019-scAPA/data/run_simulation/simulation_fq/simulation_fq_no_noise/sierra"
# )

source('/home/zhouran/data/proj/2021-0918-apa_evaluation/bin/count_peak.R')

library(Sierra)
args <- commandArgs(T)

bamfile <- normalizePath(args[1])
gtf <- normalizePath(args[2])
workdir <- args[3]
setwd(workdir)

system(command = glue::glue('ln -s {bamfile} merged.bam'))
system(command = glue::glue('ln -s {bamfile}.bai merged.bam.bai'))

regtools <- '/home/zhouran/data/proj/2021-0918-apa_evaluation/02.evaluation/Sierra/code/regtools-master/build/regtools'
reference.file <- gtf
system( command = glue::glue('{regtools} junctions extract -s 1 {bamfile} -o {bamfile}.bed'),
  wait=T)

junctions.file <-
  glue::glue("{bamfile}.bed")

# whitelist.bc.file <-
#   c(
#     paste0(extdata_path, "/example_TIP_sham_whitelist_barcodes.tsv"),
#     paste0(extdata_path, "/example_TIP_MI_whitelist_barcodes.tsv")
#   )

library(Sierra)
peak.output.file <- c("simu.txt")

FindPeaks(output.file = peak.output.file[1],   # output filename
          gtf.file = reference.file,           # gene model as a GTF file
          bamfile = bamfile[1],                # BAM alignment filename.
          junctions.file = junctions.file,     # BED filename of splice junctions exising in BAM file. 
          ncores = 4)                          # number of cores to use




peak.dataset.table = data.frame(
  Peak_file = peak.output.file,
  Identifier = "simu",
  stringsAsFactors = FALSE
)

peak.merge.output.file <- peak.output.file
count.dirs <- c("simu")
write.table('AAA', file ='barcode.tsv', quote=F,row.names=F,col.names=F)

#sham data set
CountPeaks_fake(peak.sites.file = peak.merge.output.file, 
           gtf.file = reference.file,
           bamfile = bamfile[1], 
           whitelist.file = 'barcode.tsv',
           output.dir = 'simu', 
           countUMI = FALSE, 
           ncores = 4)

## load sierra results
require(Matrix)
mtx <- Matrix::readMM('simu/matrix.mtx.gz')
pa <- read.table('simu/sitenames.tsv.gz', stringsAsFactor = F)
rownames(mtx) <- pa$V1
mtx_tmp <- data.frame(as.matrix(mtx))
colnames(mtx_tmp) <- 'counts'
saveRDS(mtx_tmp, file = 'sierra.Rds')



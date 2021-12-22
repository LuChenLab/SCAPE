args <- commandArgs(T)
bamfile = normalizePath(args[1])
workdir = args[2]
setwd(workdir)

# bamfile = '/home/zhouran/data/proj/2021-0918-apa_evaluation/03.simulation_data_with_polyA/align/simu.Aligned.sortedByCoord.out.bam'

system(command = glue::glue('ln -s {bamfile} merged.bam'))
system(command = glue::glue('ln -s {bamfile}.bai merged.bam.bai'))

# for peaks calling
# samtools view -f 128 -hb ../simu/simu.Aligned.sortedByCoord.out.bam -o simu_r2.bam
# samtools index simu_r2.bam

system(command = 'echo \'/home/zhouran/data/anaconda3/envs/apa_simu/bin/makeTagDirectory Tagdirectory <(samtools view -h merged.bam)\'|bash',wait = T)

# "findPeaks ./temp/Tagdirectory -size 50 -fragLength 100 -minDist 1 -strand separate -o ",

system(command = '/home/zhouran/data/anaconda3/envs/apa_simu/bin/findPeaks ./Tagdirectory -size 50 -fragLength 150 -minDist 1 -strand separate -o Peakfile',wait = T)


# R code
require(package = "dplyr", warn.conflicts = F)
require(package = "tidyr", warn.conflicts = F)
require(package = "ggplot2", warn.conflicts = F)
require(package = "EnvStats", warn.conflicts = F)
require(package = "parallel", warn.conflicts = F)
require(package = "mclust", warn.conflicts = F)
require("Rsubread")
require("scAPA")
bedtools.path <- '/home/zhouran/data/soft/bedtools/current/'

merge_peaks(bedtools.path = '/home/zhouran/data/soft/bedtools/current/', peaks.file = "Peakfile", path = "./")
peaks.bed <- intersect_peaks(org = 'Mm', bed.name = "./merge.peakfile.bed",
                             path = "", bedtools.path = '/home/zhouran/data/soft/bedtools/current/',
                             introns = F)
write.bed(.x = peaks.bed, f = "./peaks.bed")


wig.plus.command <- paste0(bedtools.path,
                           "bedtools genomecov -ibam merged.bam ",
                           "-bg -strand + | awk 'BEGIN ",
                           "{OFS = \"\t\"}{print $1",
                           ", $2, $3, $4, \".\", \"+\"}' > merged.wig")
system(command = wig.plus.command, wait = T)

wig.minus.command <- paste0(bedtools.path, "bedtools genomecov -ibam merged.bam ",
                            "-bg -strand - | awk 'BEGIN {OFS = \"\t\"}{print $1",
                            ", $2, $3, $4, \".\", \"-\"}' >> merged.wig")
system(command = wig.minus.command, wait = T)

intersect.wig.command <- paste0(bedtools.path, "bedtools intersect -s -wb ",
                                "-b peaks.bed -a merged.wig > intersected.wig")
system(intersect.wig.command, wait = T)
peaks.wig <- read.delim(file = "intersected.wig", header = F)
peaks.wig <- split(x = peaks.wig, f = peaks.wig$V10, drop = T)


mclust.command <- paste0("Mclust()")
cores <- 10

bed <- plyr::rbind.fill(parallel::mclapply(1:length(peaks.wig),
                                           FUN = creat_mclus,
                                           mc.cores = cores,
                                           mc.preschedule = T))


utr.saf <- bed[, c(4, 1, 2, 3, 6)]
rm("bed")
colnames(utr.saf) <- c("GeneID", "Chr", "Start", "End", "Strand")



clusters <- 'merged'
bam.cluster.files <- paste0(clusters, ".bam")

counts <- Rsubread::featureCounts(files = bam.cluster.files, isGTFAnnotationFile = F,
                        strandSpecific = 1, annot.ext = utr.saf,
                        largestOverlap = T, nthreads = cores)


co <- cbind.data.frame(rownames(counts$counts), counts$counts)
colnames(co) <- c("Peak_ID", clusters)
meta <- counts$annotation
meta <- meta[, c(2, 3, 4, 1, 6, 5)]

final_res <- merge(co, utr.saf,by.x='Peak_ID', by.y='GeneID')
final_res$ensembl <- sapply(strsplit(as.character(final_res$Peak_ID), split='[.]'),'[[',1)
library(org.Mm.eg.db)
# bitr(geneID, fromType, toType, OrgDb, drop = TRUE)

id_map <- clusterProfiler::bitr(unique(final_res$ensembl), fromType='ENSEMBL', toType='SYMBOL',OrgDb=org.Mm.eg.db)
# truth <- readRDS('../groud_t.Rds')
# id_map <- id_map[id_map$SYMBOL %in% names(truth),]

final_res$gene_id <- NA
final_res$gene_id <- plyr::mapvalues(from = id_map$ENSEMBL, to = id_map$SYMBOL, x=final_res$ensembl)

# 这里还要根据Peak_ID来获得3‘的最后位点，来确定最后的结果。
final_res_collapse <- lapply(split(final_res,final_res$Peak_ID), function(x){
  if (dim(x)[1] == 1){
    return(x)
  }
  if (x$Strand[1] =='+') {
    return(x[x$End==max(x$End),])
  } else {
    return(x[x$Start==min(x$Start),])
  }
})

final_res_collapse <- do.call(rbind, final_res_collapse)
saveRDS(final_res_collapse, file = 'scAPA.raw.Rds')

# fasta.path <- '/mnt/raid61/Microwell/mm10/fasta/genome.fa'
# char.length.path <- '/mnt/raid61/Microwell/mm10/fasta/genome.fa.fai'
# metadata <- read_down.seq(saf = utr.saf,
#                           char.length.path = char.length.path,
#                           fasta.path = fasta.path, chr.modify = F)

# aseq <- metadata[, c(4, 6)]
# a <- set_scAPAList(.clus.counts = co, .row.Data = meta, .down.seq = aseq)
# saveRDS(object = a, file = "../outs/Peaks.RDS")
# write_log(f = "../scAPA.script.log", stage = "featureCounts for 3UTRs")
# if (int) {
#     write_log_start(f = "../scAPA.script.log",
#                     "featureCounts for introns", command = NA)
#     counts_int <- Rsubread::featureCounts(files = bam.cluster.files,
#                                 isGTFAnnotationFile = F,
#                                 annot.ext = int.saf, strandSpecific = 1,
#                                 largestOverlap = T, nthreads = c)
#     co_int <- cbind.data.frame(rownames(counts_int$counts),
#                                counts_int$counts)
#     colnames(co_int) <- c("Peak_ID", clusters)
#     meta_int <- counts_int$annotation
#     meta_int <- meta_int[, c(2, 3, 4, 1, 6, 5)]
#     metadata_int <- read_down.seq(saf = int.saf,
#                                   char.length.path = char.length.path,
#                                   fasta.path = fasta.path, chr.modify = T)
#     aseq_int <- metadata_int[, c(4, 6)]
#     a.int <- set_scAPAList(.clus.counts = co_int, .row.Data = meta_int,
#                            .down.seq = aseq_int)
#     saveRDS(object = a.int, file = "../outs/Peaks.int.RDS")
#     write_log(f = "../scAPA.script.log",
#               stage = "featureCounts for introns")
# }


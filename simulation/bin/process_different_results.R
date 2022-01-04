rm(list = ls())
library(data.table)
library(tibble)
library(dplyr)
core_use <- 12

args <- commandArgs(T)

collapse_file <- args[1]
truth <- args[2]
##################### ground truth results
truth <-
  read.table(
    truth,
    stringsAsFactors = F
  )

truth$truth <- truth$V2 * truth$V3
colnames(truth) <-
  c('meta', 'coverage', 'pa_ws', 'noise_ws', 'pa_coverage')

collapse_info <-
  read.table(
    collapse_file,
    stringsAsFactors = F
  )

collapse_info$region <- paste(collapse_info$V1,
                              collapse_info$V2,
                              collapse_info$V3,
                              collapse_info$V4,
                              sep = '_')

collapse_info$pasite <-
  ifelse(collapse_info$V4 == "+", collapse_info$V7, collapse_info$V6)

collapse_info <- collapse_info[, c('V8', 'region', 'pasite')]
colnames(collapse_info)[1] <- 'meta'


truth <- merge(truth, collapse_info, by = 'meta')

truth <- split(truth, truth$region)

truth <- lapply(truth, function(x) {
  # x <- truth[[3]]
  if (isTRUE(dim(x)[1] == 1)) {
    return(x)
  } else {
    x$new_ws <- x$coverage / sum(x$coverage)
    return(x)
  }
})

gene_id <- unname(unlist(lapply(truth, function(x) {
  sapply(strsplit(x$meta[[1]], split = '_|,'), '[[', 2)
})))
names(truth) <- gene_id

saveRDS(truth, file = 'ground_truth.Rds')


##################### MAPPER results
dat <- readRDS('MAAPER/result.rds')

parallel::mclapply(names(dat), function(x) {
  if (is.null(dat[[x]])) {
    return(NULL)
  }
  # message(x)
  
  if ('pas' %in% names(dat[[x]])) {
    pa_site <- gsub('^chr', '', dat[[x]]$pas)
    pa_weight <- 1
    pa_counts <- pa_weight * dat[[x]]$nreads$c1
  } else{
    pa_site <- gsub('^chr', '', rownames(dat[[x]]$alp_alt))
    pa_weight <- dat[[x]]$alp_alt[, 1]
    pa_counts <- pa_weight * dat[[x]]$nreads$c1
  }
  
  tibble(
    gene_id = x,
    chrom = sapply(strsplit(pa_site, split = ':'), '[[', 1),
    start = sapply(strsplit(pa_site, split = ':'), '[[', 2),
    end = sapply(strsplit(pa_site, split = ':'), '[[', 2),
    strand = sapply(strsplit(pa_site, split = ':'), '[[', 3),
    pa_site = pa_site,
    pa_weight = pa_weight,
    pa_counts = pa_counts,
    pa_id = paste(gene_id, chrom, start, end, ifelse(strand == '+', 1,-1), sep = ':'),
    pa = sapply(strsplit(pa_site, split = ':'), '[[', 2)
  )
}, mc.cores = core_use) -> dat

dat <- do.call(rbind, dat)

saveRDS(dat, file = 'MAAPER/MAAPER.Rds')



##################### SCAPTURE results
# SCAPTURE's results like sierra
rm(list = ls())
dat <- data.table::fread('SCAPTURE/simu.KeepCell.assigned')

dat$Chr <-
  lapply(strsplit(dat$Chr, split = ';'), unique) %>% unlist()
dat$Strand <-
  lapply(strsplit(dat$Strand, split = ';'), unique) %>% unlist()

dat$Start <-
  lapply(strsplit(dat$Start, split = ';'), function(x) {
    min(as.numeric(x))
  }) %>% unlist()

dat$End <-
  lapply(strsplit(dat$End, split = ';'), function(x) {
    max(as.numeric(x))
  }) %>% unlist()

colnames(dat)[7] <- 'counts'
dat$gene_name <- sapply(strsplit(dat$Geneid, split = '-'), '[[', 1)
# Wars2:3:99218517-99218865:1
dat$paid <- paste(dat$gene_name,
                  dat$Chr,
                  dat$Start,
                  dat$End,
                  ifelse(dat$Strand == '+', 1, -1),
                  sep = ':')
dat$pa_site <-
  paste(dat$Chr,
        ifelse(dat$Strand == '+', dat$End, dat$Start),
        dat$Strand,
        sep = ':')


pa_weights <- dat %>% group_by(gene_name) %>%
  mutate(pa_weights = counts / sum(counts)) %>%
  ungroup() %>%
  pull(pa_weights)

dat <- data.table(
  gene_id = dat$gene_name,
  chrom = dat$Chr,
  start = dat$Start,
  end = dat$End,
  strand = dat$Strand,
  pa_site = dat$pa_site,
  pa_weights = pa_weights,
  pa_counts = dat$counts,
  pa_id = dat$paid,
  pa = ifelse(dat$Strand == '+', dat$End, dat$Start)
)

saveRDS(dat, file = 'SCAPTURE/SCAPTURE.Rds')

##################### scAPA results
rm(list = ls())
dat <- readRDS("scAPA/scAPA.raw.Rds")

dat$pa_site <-
  paste(dat$Chr,
        ifelse(dat$Strand == '+', dat$End, dat$Start),
        dat$Strand,
        sep = ':')
dat$paid <- paste(dat$gene_id,
                  dat$Chr,
                  dat$Start,
                  dat$End,
                  ifelse(dat$Strand == '+', 1,-1),
                  sep = ':')
pa_weights <- dat %>% group_by(gene_id) %>%
  mutate(pa_weights = merged / sum(merged)) %>%
  ungroup() %>%
  pull(pa_weights)

dat <- data.table(
  gene_id = dat$gene_id,
  chrom = dat$Chr,
  start = dat$Start,
  end = dat$End,
  strand = dat$Strand,
  pa_site = dat$pa_site,
  pa_weights = pa_weights,
  pa_counts = dat$counts,
  pa_id = dat$paid,
  pa = ifelse(dat$Strand == '+', dat$End, dat$Start)
)

saveRDS(dat, file = 'scAPA/scAPA.Rds')

##################### scAPAtrap results
rm(list = ls())
library(dplyr)
library(data.table)
library(ChIPseeker)
library(GenomicFeatures)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
txdbmm10 <- movAPA::parseGenomeAnnotation(txdb)
trap <- fread('scAPAtrap/peak_assigned')
trap$coord <- ifelse(trap$Strand == '+', trap$End, trap$Start)
trap$id <- paste(trap$Chr, trap$coord, trap$Strand, sep = ':')

tmp <- trap[, c('Chr', 'coord', 'Strand')]
colnames(tmp) <- c('chr', 'coord', 'strand')
tmp <- movAPA::annotatePAC(tmp, aGFF = txdbmm10)

test <- clusterProfiler::bitr(tmp$gene,
                              fromType = 'ENTREZID',
                              toType = 'SYMBOL',
                              OrgDb = org.Mm.eg.db)
tmp$gene_id <-
  plyr::mapvalues(from = test$ENTREZID,
                  to = test$SYMBOL,
                  x = tmp$gene)
tmp <- tmp[tmp$gene_id != "character(0)"]

tmp$id <-
  paste(tmp$chr, tmp$coord, tmp$strand, sep = ':')

trap$gene_id <-
  plyr::mapvalues(
    from = tmp$id,
    to = tmp$gene_id,
    x = trap$id,
    warn_missing = F
  )

trap$gene_id[trap$id == trap$gene_id] <- 'Unkown_gene'

pa_weights <- trap %>% group_by(gene_id) %>%
  mutate(pa_weights = merged.bam / sum(merged.bam)) %>%
  ungroup() %>%
  pull(pa_weights)

trap$paid <- paste(trap$gene_id,
                   trap$Chr,
                   trap$Start,
                   trap$End,
                   ifelse(trap$Strand == '+', 1,-1),
                   sep = ':')

trap$pa_site <-
  paste(trap$Chr,
        trap$coord,
        trap$Strand,
        sep = ':')

dat <- data.table(
  gene_id = trap$gene_id,
  chrom = trap$Chr,
  start = trap$Start,
  end = trap$End,
  strand = trap$Strand,
  pa_site = trap$pa_site,
  pa_weights = pa_weights,
  pa_counts = trap$merged.bam,
  pa_id = trap$paid,
  pa = trap$coord
)

dat$chrom <- gsub('^chr', '', dat$chrom)
dat$pa_site <- gsub('^chr', '', dat$pa_site)
dat$pa_id <- gsub(':chr', ':', dat$pa_id)

saveRDS(dat, file = 'scAPAtrap/scAPAtrap.Rds')

##################### Sierra results
rm(list = ls())
dat <- readRDS('Sierra/sierra.Rds')

dat <- cbind(dat, data.table(do.call(rbind, strsplit(
  rownames(dat), split = ':'
))))

dat$gene_name <- dat$V1
dat$Chr <- dat$V2
dat$Strand <- ifelse(dat$V4 == 1 , '+', '-')
dat$Start <- sapply(strsplit(dat$V3, split = '-'), '[[', 1)
dat$End <- sapply(strsplit(dat$V3, split = '-'), '[[', 2)
dat$paid <- paste(dat$gene_name,
                  dat$Chr,
                  dat$Start,
                  dat$End,
                  ifelse(dat$Strand == '+', 1,-1),
                  sep = ':')


pa_weights <- dat %>% group_by(gene_name) %>%
  mutate(pa_weights = counts / sum(counts)) %>%
  ungroup() %>%
  pull(pa_weights)
dat$pa_site <-
  paste(dat$Chr,
        ifelse(dat$Strand == '+', dat$End, dat$Start),
        dat$Strand,
        sep = ':')


dat <- data.table(
  gene_id = dat$gene_name,
  chrom = dat$Chr,
  start = dat$Start,
  end = dat$End,
  strand = dat$Strand,
  pa_site = dat$pa_site,
  pa_weights = pa_weights,
  pa_counts = dat$counts,
  pa_id = dat$paid,
  pa = ifelse(dat$Strand == '+', dat$End, dat$Start)
)

saveRDS(dat, file = 'Sierra/Sierra.Rds')




##################### scDapars results

# dat <- data.table(
#   gene_id = dat$gene_name,
#   chrom = dat$Chr,
#   start = dat$Start,
#   end = dat$End,
#   strand = dat$Strand,
#   pa_site = dat$pa_site,
#   pa_weights = pa_weights,
#   pa_counts = dat$counts,
#   pa_id = dat$paid
# )


dapars2 <- read.table('scDapars/dapars2_res.txt',
                      stringsAsFactors = F,
                      header = T)
dapars2$Strand <-
  sapply(strsplit(unlist(sapply(
    strsplit(dapars2$Gene, split = ','), tail, 1
  )), split = '_'), '[[', 7)

# dapars2$Loci <- gsub('-', '_', dapars2$Loci)
# dapars2$region <-
#   gsub(':', '_', paste(dapars2$Loci, dapars2_strand, sep = '_'))
dapars2$gene_id <-
  sapply(strsplit(dapars2$Gene, split = ',|_'), '[[', 2)
dapars2$Chr <-
  sapply(strsplit(dapars2$Loci, split = ':|-'), '[[', 1)

dapars2$Start <-
  dapars2$Predicted_Proximal_APA
dapars2$End <-
  dapars2$Predicted_Proximal_APA

dapars2$pa_site <-
  paste(dapars2$Chr,
        dapars2$Predicted_Proximal_APA,
        dapars2$Strand,
        sep = ':')
dapars2$paid <- paste(
  dapars2$gene_id,
  dapars2$Chr,
  dapars2$Start,
  dapars2$End,
  ifelse(dapars2$Strand == '+', 1, -1),
  sep = ':'
)

dat <- data.table(
  gene_id = dapars2$gene_id,
  chrom = dapars2$Chr,
  start = dapars2$Start,
  end = dapars2$End,
  strand = dapars2$Strand,
  pa_site = dapars2$pa_site,
  pa_weights = dapars2[,grep('PDUI', colnames(dapars2))],
  pa_counts = NA,
  pa_id = dapars2$paid,
  pa = dapars2$Predicted_Proximal_APA
)

saveRDS(dat, file = 'scDapars/scDapars.Rds')

##################### polyApipe results
rm(list = ls())

dat <- data.table::fread('polyApipe/simu_counts.tab.gz')
dat_to_annot <-
  data.table(do.call(rbind, strsplit(dat$gene, split = '_')))
dat_to_annot$V3 <- ifelse(dat_to_annot$V3 == 'f', '+', '-')
dat_to_annot$V1 <- paste0('chr', dat_to_annot$V1)
colnames(dat_to_annot) <- c('chr', 'coord', 'strand')

dat$id <-
  paste(dat_to_annot$chr,
        dat_to_annot$coord,
        dat_to_annot$strand,
        sep = ':')
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
txdbmm10 <- movAPA::parseGenomeAnnotation(txdb)

tmp <- movAPA::annotatePAC(dat_to_annot, aGFF = txdbmm10)

test <- clusterProfiler::bitr(tmp$gene,
                              fromType = 'ENTREZID',
                              toType = 'SYMBOL',
                              OrgDb = org.Mm.eg.db)
tmp$gene_id <-
  plyr::mapvalues(from = test$ENTREZID,
                  to = test$SYMBOL,
                  x = tmp$gene)
tmp <- tmp[tmp$gene_id != "character(0)"]

tmp$id <-
  paste(tmp$chr, tmp$coord, tmp$strand, sep = ':')


dat$gene_id <-
  plyr::mapvalues(
    from = tmp$id,
    to = tmp$gene_id,
    x = dat$id,
    warn_missing = F
  )

dat$gene_id[dat$id == dat$gene_id] <- 'Unkown_gene'

pa_weights <- dat %>% group_by(gene_id) %>%
  mutate(pa_weights = count / sum(count)) %>%
  ungroup() %>%
  pull(pa_weights)

dat$paid <- paste(
  dat$gene_id,
  sapply(strsplit(dat$gene, split = '_'), '[[', 1),
  sapply(strsplit(dat$gene, split = '_'), '[[', 2),
  sapply(strsplit(dat$gene, split = '_'), '[[', 2),
  ifelse(sapply(strsplit(dat$gene, split = '_'), '[[', 2) == 'f', 1, -1),
  sep = ':'
)

dat <- data.table(
  gene_id = dat$gene_id,
  chrom = sapply(strsplit(dat$gene, split = '_'), '[[', 1),
  start = sapply(strsplit(dat$gene, split = '_'), '[[', 2),
  end = sapply(strsplit(dat$gene, split = '_'), '[[', 2),
  strand = sapply(strsplit(dat$id, split = ':'), '[[', 3),
  pa_site = gsub('^chr', '', dat$id),
  pa_weights = pa_weights,
  pa_counts = dat$count,
  pa_id = dat$paid,
  pa = sapply(strsplit(dat$gene, split = '_'), '[[', 2)
)

saveRDS(dat, file = 'polyApipe/polyApipe.Rds')


## SCAPE

rm(list = ls())

dat <- data.table::fread('pasite.csv.gz')
colnames(dat)[2] <- 'count'
dat_to_annot <-
  data.table(do.call(rbind, strsplit(dat$V1, split = ':')))
dat_to_annot$V1 <- paste0('chr', dat_to_annot$V1)
dat_to_annot <- dat_to_annot[, c(1, 2, 4)]
colnames(dat_to_annot) <- c('chr', 'coord', 'strand')

dat$id <-
  paste(dat_to_annot$chr,
        dat_to_annot$coord,
        dat_to_annot$strand,
        sep = ':')
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
txdbmm10 <- movAPA::parseGenomeAnnotation(txdb)

tmp <- movAPA::annotatePAC(dat_to_annot, aGFF = txdbmm10)

test <- clusterProfiler::bitr(tmp$gene,
                              fromType = 'ENTREZID',
                              toType = 'SYMBOL',
                              OrgDb = org.Mm.eg.db)
tmp$gene_id <-
  plyr::mapvalues(from = test$ENTREZID,
                  to = test$SYMBOL,
                  x = tmp$gene)
tmp <- tmp[tmp$gene_id != "character(0)"]

tmp$id <-
  paste(tmp$chr, tmp$coord, tmp$strand, sep = ':')


dat$gene_id <-
  plyr::mapvalues(
    from = tmp$id,
    to = tmp$gene_id,
    x = dat$id,
    warn_missing = F
  )

dat$gene_id[dat$id == dat$gene_id] <- 'Unkown_gene'



pa_weights <- dat %>% group_by(gene_id) %>%
  mutate(pa_weights = count / sum(count)) %>%
  ungroup() %>%
  pull(pa_weights)

dat$paid <- paste(
  dat$gene_id,
  sapply(strsplit(dat$V1, split = ':'), '[[', 1),
  sapply(strsplit(dat$V1, split = ':'), '[[', 2),
  sapply(strsplit(dat$V1, split = ':'), '[[', 2),
  sapply(strsplit(dat$V1, split = ':'), '[[', 4),
  sep = ':'
)

dat <- data.table(
  gene_id = dat$gene_id,
  chrom = sapply(strsplit(dat$V1, split = ':'), '[[', 1),
  start = sapply(strsplit(dat$V1, split = ':'), '[[', 2),
  end = sapply(strsplit(dat$V1, split = ':'), '[[', 2),
  strand = sapply(strsplit(dat$V1, split = ':'), '[[', 4),
  pa_site = gsub('^chr', '', dat$id),
  pa_weights = pa_weights,
  pa_counts = dat$count,
  pa_id = dat$paid,
  pa = sapply(strsplit(dat$V1, split = ':'), '[[', 2)
)

saveRDS(dat, file = 'SCAPE/SCAPE.Rds')

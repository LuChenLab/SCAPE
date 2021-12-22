library(scAPAtrap)
args <- commandArgs(T)
bamfile <- normalizePath(args[1])
workdir <- args[2]
setwd(workdir)

thread <- 4
samtools.path <- "samtools"
featureCounts.path <- "/home/zhouran/data/soft/featureCounts/current/bin/featureCounts"


system(command = glue::glue('ln -s {bamfile} merged.bam'))
system(command = glue::glue('ln -s {bamfile}.bai merged.bam.bai'))

input <- "merged.bam"
outputF <- paste0(gsub(".bam", "", input), ".forward.bam")
outputR <- paste0(gsub(".bam", "", input), ".reverse.bam")

samtools.Fcommand <- paste0(samtools.path, " view -@ ", thread,
    " -h -F 0x10 -bS ", input, " > ", outputF)
system(command = samtools.Fcommand, wait = T)
samtools.Rcommand <- paste0(samtools.path, " view -@ ", thread,
    " -h -f 0x10 -bS ", input, " > ", outputR)
system(command = samtools.Rcommand, wait = T)

command <- paste0(samtools.path, " index -@ ", thread, " ",
    outputF)
system(command = command, wait = T)
command <- paste0(samtools.path, " index -@ ", thread, " ",
    outputR)
system(command = command, wait = T)

nextinput <- c('merged.forward.bam', 'merged.reverse.bam')
chrs <- c(as.character(1:19),'X','Y')
maxwidth <- 1000
readlength <- 98
outputdir <- './result'
fullcovF <- loadBpCoverages(nextinput[1],chrs)
fullcovR <- loadBpCoverages(nextinput[2],chrs)

forwardPeaks <-findPeaks(fullcovF, '+', readlength, maxwidth)
reversePeaks <-findPeaks(fullcovR, '-', readlength, maxwidth)

head(forwardPeaks)
head(reversePeaks)

peaksfile <- generateSAF(forwardPeaks, reversePeaks, outputdir)
peaksfile

# here to generate final `peak_assigned` file.
final.bam <- generateFinalBam(featureCounts.path,samtools.path,input,peaksfile,24)

# final.bam <- 'final.bam'
# umitools.path <- '/home/zhouran/data/anaconda3/envs/polyapipe/bin/umi_tools'
# counts.tsv <- countPeaks(umitools.path,final.bam,outputdir,TenX=T)
# counts.tsv <- 'result/counts.tsv.gz'

findChrTails_2 <- function(bamfile, chr, len) {
  which <- GenomicRanges::GRanges(data.frame(
    seqnames = chr,
    start = 1,
    end = len
  ))
  
  what <- c("rname", "strand", "pos", "cigar", "seq")
  param <- Rsamtools::ScanBamParam(what = what, which = which)
  gal1 <-
    GenomicAlignments::readGAlignments(bamfile, use.names = TRUE, param = param)
  s_1 <-
    (grepl("[0-9]*M[1-9]{2,}S", gal1@cigar) &
       as.vector(gal1@strand) == "+")
  s_2 <-
    (grepl("[1-9]{2,}S[0-9]*M", gal1@cigar) &
       as.vector(gal1@strand) == "-")
  
  bam1 <- gal1[s_1]
  bam2 <- gal1[s_2]
  
  bam1 <-
    bam1[grepl(
      "(A{3,}[^A]{0,2})*A{6,}([^A]{0,2}A{3,})*.{0,2}?$",
      bam1@elementMetadata@listData$seq
    )]
  bam2 <-
    bam2[grepl(
      "^.{0,2}?(T{3,}[^T]{0,2})*T{6,}([^T]{0,2}T{3,})*",
      bam2@elementMetadata@listData$seq
    )]
  
  final_bam1 <-
    data.frame(
      chr = as.vector(seqnames(bam1)),
      strand = as.vector(strand(bam1)),
      coord = end(bam1)
    )
  final_bam2 <-
    data.frame(
      chr = as.vector(seqnames(bam2)),
      strand = as.vector(strand(bam2)),
      coord = start(bam2)
    )
  
  bam <- rbind(final_bam1, final_bam2)
  bam <-
    dplyr::group_by(bam, chr, strand, coord) %>% summarise(count = n())
  
  return(bam)
}

abamfile <- 'merged.bam'
samtoolsexc <- '/usr/bin/samtools'
chrID<-system( paste(samtoolsexc,"view -H", abamfile ," | grep '^@SQ' | cut -f2 "), intern = T)
chrID<-lapply(chrID, function(x){
  x<-strsplit(x, ":")[[1]][[2]]
})

chrL<-system( paste(samtoolsexc,"view -H", abamfile ," | grep '^@SQ' | cut -f3 "), intern = T)
chrL<-lapply(chrL, function(x){
  x<-as.integer(strsplit(x, ":")[[1]][[2]])
})


lapply(1:length(chrID), function(x){
    tails <- findChrTails_2(bamfile = 'merged.bam', chr=chrID[[x]], len = chrL[[x]])   
}) -> tails

tails <- do.call(rbind, tails)
saveRDS(tails, file = 'tails.Rds')


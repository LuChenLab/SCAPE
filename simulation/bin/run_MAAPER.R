args <- commandArgs(T)

pas_annotation <- args[1]
gtf <- args[2]
outdir <- args[3]
bam_c1 <- args[4]
library(MAAPER)

pas_annotation = readRDS(pas_annotation)
# gtf = "/mnt/raid/Ref/MusMus/release101/mm10/genes/genes.gtf"
# # bam file of condition 1 (could be a vector if there are multiple samples)
# bam_c1 = "/home/zhouran/data/proj/2021-0918-apa_evaluation/03.simulation_data_with_polyA/align/simu.Aligned.sortedByCoord.out.bam"
# # bam file of condition 2 (could be a vector if there are multiple samples)

maaper(gtf, 
       pas_annotation, 
       output_dir = outdir, 
       bam_c1, bam_c1,
       read_len = 150, 
       ncores = 12  
      )

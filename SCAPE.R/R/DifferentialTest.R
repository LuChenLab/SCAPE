
#'@title Find differentail expression APA events
#'@description Find DE APA events by DEXseq test
#'@param obj A Seurat object which contained apa matrix information.
#'@param idents.1 A sub-cluster from Idents(dat)
#'@param idents.2 A sub-cluster from Idents(dat). If idents.2 == NULL, then idents.1 versus other all cells.
#'@param annot The annotation of all pA sites, generate from AnnotationSite.
#'@param paLoc The genomic annotation of pA for performing DE test.
#'@param slot The slot data of Seurat, including counts, data and scale.data. defualt: counts
#'@param assay The name of assay in Seurat class, default: RNA
#'@param seedUse The random seed used.
#'@param cores The num of cpu for DE test.
#'@param num.split The num of cell chunks for random sample from all cells.
#'#'@example
#'\dontrun{test <- FindDE(
#'dat,
#'idents.1 = levels(Idents(dat))[1],
#'idents.2 = levels(Idents(dat))[2],
#'annot = pa_anno,
#'assay = "RNA",
#'cores = 20
#')
#'}
#'@export


FindDE <- function(obj,
                   idents.1,
                   idents.2=NULL,
                   annot,
                   paLoc = '3UTRs',
                   slot = 'counts',
                   assay = 'RNA',
                   seedUse = 1,
                   cores = 1,
                   num.splits = 6) {
  set.seed(seedUse)
  # here to keep these pa which located in 3UTRs, and was part of APA evens.
  annot <-
    annot[annot$pa_site %in% rownames(GetAssayData(obj, slot = 'counts', assay = assay)) &
            !is.na(annot$annot.symbol) &
            annot$Type %in% paLoc, ]

  freq_to_filter <- table(annot$annot.symbol)

  annot <-
    annot[annot$annot.symbol %in% names(freq_to_filter)[freq_to_filter != 1],]

  pa_use <- annot$pa_site

  cells.1 = names(Seurat::Idents(obj))[which(Seurat::Idents(obj) == idents.1)]
  if (is.null(idents.2)) {
    cells.2 = names(Seurat::Idents(obj))[which(Seurat::Idents(obj) != idents.1)]
  } else {
    cells.2 = names(Seurat::Idents(obj))[which(Seurat::Idents(obj) == idents.2)]
  }


  cells.sub.1 = split(cells.1, sort(1:length(cells.1) %% num.splits))
  data_to_test.1 = matrix(, nrow = length(pa_use), ncol = length(cells.sub.1))

  cells.sub.2 = split(cells.2, sort(1:length(cells.2) %% num.splits))
  data_to_test.2 = matrix(, nrow = length(pa_use), ncol = length(cells.sub.2))


  for (i in 1:length(cells.sub.1)) {
    this.set <- cells.sub.1[[i]]
    sub.matrix <-
      GetAssayData(obj, slot = "counts", assay = assay)[pa_use, this.set]
    if (length(this.set) > 1) {
      this.profile <- as.numeric(apply(sub.matrix, 1, function(x)
        sum(x)))
      data_to_test.1[, i] <- this.profile
    } else {
      data_to_test.1[, i] <- sub.matrix
    }
  }

  rownames(data_to_test.1) <- pa_use
  colnames(data_to_test.1) <-
    paste0("Population1_", 1:length(cells.sub.1))

  for (i in 1:length(cells.sub.2)) {
    this.set <- cells.sub.2[[i]]
    sub.matrix <-
      Seurat::GetAssayData(obj, slot = "counts", assay = assay)[pa_use, this.set]
    if (length(this.set) > 1) {
      this.profile <- as.numeric(apply(sub.matrix, 1, function(x)
        sum(x)))
      data_to_test.2[, i] <- this.profile
    } else {
      data_to_test.2[, i] <- sub.matrix
    }
  }

  rownames(data_to_test.2) <- pa_use
  colnames(data_to_test.2) <-
    paste0("Population2_", 1:length(cells.sub.1))


  peak.matrix <- cbind(data_to_test.1, data_to_test.2)
  rownames(peak.matrix) <-
    stringr::str_replace_all(rownames(peak.matrix), ':', '_')

  sampleTable <-
    data.frame(row.names = c(colnames(data_to_test.1), colnames(data_to_test.2)),
               condition = c(rep("target", ncol(data_to_test.1)),
                             rep("comparison", ncol(data_to_test.2))))

  exon_info <- annot[, c('pa_site', 'annot.symbol')]
  exon_info$pa_site <-
    stringr::str_replace_all(exon_info$pa_site, ':', '_')


  dxd <- DEXSeq::DEXSeqDataSet(
    peak.matrix[, rownames(sampleTable)],
    sampleTable,
    design = ~ sample + exon + condition:exon,
    exon_info$pa_site,
    exon_info$annot.symbol,
    featureRanges = NULL,
    transcripts = NULL,
    alternativeCountData = NULL
  )
  require(magrittr)
  if (cores != 1) {
    BPPARAM = BiocParallel::MulticoreParam(cores)
    dxd %<>%
      estimateSizeFactors %<>%
      estimateDispersions(BPPARAM = BPPARAM) %<>%
      testForDEU(BPPARAM = BPPARAM) %<>%
      estimateExonFoldChanges(BPPARAM = BPPARAM)

  } else {
    dxd %<>%
      estimateSizeFactors %<>%
      estimateDispersions %<>%
      testForDEU %<>%
      estimateExonFoldChanges
  }
  dxr1 = DEXSeqResults(dxd)
  dxr1 <-
    dxr1[, c("groupID",
             "pvalue",
             "padj",
             "log2fold_target_comparison")]
  dxr1 <- as.data.frame(dxr1)
  pa_use <-
    stringr::str_replace_all(sapply(strsplit(rownames(dxr1), split = ":"), "[[", 2),
                             "_", ":")
  colnames(dxr1) <- c('gene', 'p_val', 'p_val_adj', 'avg_logFC')
  dxr1$pa_site <- pa_use

  dxr1[["pct.1"]] <-
    apply(
      X = GetAssayData(obj, "counts", assay = assay)[pa_use, Idents(obj) == idents.1],
      MARGIN = 1,
      FUN = Seurat:::PercentAbove,
      threshold = 0
    )

  if (is.null(idents.2)) {
    dxr1[["pct.2"]] <-
      apply(
        X = GetAssayData(obj, "counts", assay = assay)[pa_use, Idents(obj) != idents.1],
        MARGIN = 1,
        FUN = Seurat:::PercentAbove,
        threshold = 0
      )
    dxr1[['versus']] <- glue::glue('{idents.1}')
  } else {
    dxr1[['versus']] <- glue::glue('{idents.1}_Vs_{idents.2}')
    dxr1[["pct.2"]] <-
      apply(
        X = GetAssayData(obj, "counts", assay = assay)[pa_use, Idents(obj) == idents.2],
        MARGIN = 1,
        FUN = Seurat:::PercentAbove,
        threshold = 0
      )

  }

  dxr1 <-
    dxr1[, c('gene',
             'pa_site',
             'pct.1',
             'pct.2',
             'versus',
             'p_val',
             'p_val_adj',
             'avg_logFC')]

  return(dxr1)
}
#'@title Differential expression test for APA events based bayes factor
#'@description Differential expression test for APA events based bayes factor
#'@param obj A Seurat object after calculated psi (psi function).
#'@param idents.1 A sub-cluster from Idents(dat)
#'@param idents.2 A sub-cluster from Idents(dat). If idents.2 == NULL, then idents.1 versus other all cells.
#'@param slot The slot data of Seurat, including counts, data and scale.data. defualt: counts
#'@param assay The name of assay in Seurat class, default: apapsi
#'@param seedUse The random seed used.
#'@param cores The num of cpu for DE test.
#'#'@example
#'\dontrun{test <- FindBF(
#'dat,
#'idents.1 = levels(Idents(dat))[1],
#'idents.2 = levels(Idents(dat))[2],
#'slot = counts,
#'assay = "apapsi",
#'cores = 20,
#'
#')
#'}
#'@export

FindBF <- function(obj,
                   idents.1,
                   idents.2 = NULL,
                   slot = 'counts',
                   assay = 'apapsi',
                   cores = 1,
                   seedUse = 1) {
  set.seed(seedUse)
  cells.1 = names(Seurat::Idents(obj))[which(Seurat::Idents(obj) == idents.1)]
  
  # here for idents.1
  mat.idents.1 <-
    GetAssayData(obj,
                 slot = slot,
                 assay = assay)[, cells.1]
  
  # here for idents.2
  if (is.null(idents.2)) {
    cells.2 = names(Seurat::Idents(obj))[which(Seurat::Idents(obj) != idents.1)]
  } else {
    cells.2 = names(Seurat::Idents(obj))[which(Seurat::Idents(obj) == idents.2)]
  }

  mat.idents.2 <- 
    GetAssayData(obj,
                 slot = slot,
                 assay = assay)[, cells.2]
  parallel::mclapply(rownames(mat.idents.1), function(x){
    val_1 <- na.omit(mat.idents.1[x, ])
    val_2 <- na.omit(mat.idents.2[x, ])
    BayesFactor::ttestBF(
      x=val_1,
      y=val_2
    )@bayesFactor$bf -> bf_val
    delta_psi <- mean(val_1) - mean(val_2)
    val <- c(bf_val, delta_psi)
    names(val)-> c('BayesFactor', 'Delta_psi')
    val
  }) -> res
  as.data.frame(do.call(rbind,res)) -> res
  res$pa_site <- rownames(mat.idents.1)
  res
}






#'@title Psi calculation.
#'@description Psi calculation.
#'@param obj A Seurat object include `apa` assay.
#'@param annot Annotation dataframe.
#'@param chunk Cell number per chunk for reducing memory usage.
#'@param cores The num of cpu for calculating psi.
#'@example
#'@example
#'@example
#'\dontrun{annot_info <-
#'   annotate_from_gtf('/home/zhouran/data/tmp/annotate_test/chr1.gtf.gz',
#'                     'Mm10',
#'                    cores = 10)
#'}
#'
#'@export

psi <- function(obj,
                annot,
                chunk = 1000,
                assay = 'apapsi',
                cores = 1) {
  # if annotation dataframe not generate from AnnotationSite.

  stopifnot("Type" %in% colnames(annot))

  mat <- Seurat::GetAssayData(obj, 'counts', assay = 'apa')
  annot <- annot %>%
    subset(Type != "INTERGENIC") %>%
    dplyr::filter(!is.na(annot.gene_id)) %>%
    dplyr::select(annot.gene_id, pa_site)
  colnames(annot) <- c('V2', 'V1')
  rownames(annot) <- annot[['V1']]

  ##keep PA.id within ensembl gene
  annoG <- annot[grep("^ENS", annot$V2), ]
  add.id <- annoG[rownames(mat),]
  add.id <- na.omit(add.id)

  pa.geneid.matrix <- mat[add.id$V1,]

  # ratio calculating
  ratio <- function(x) {
    if (is.numeric(x)) {
      n = x / sum(x)
    } else{
      n = x
    }
    return(n)
  }
  #convert NA to 0
  DT.for.set.sqln  <- function(x) {
    for (j in seq_len(ncol(x)))
      data.table::set(x, which(is.na(x[[j]])), j, 0)
  }
  PAratio <- function(mtx) {
    mtx_dt <- data.table::as.data.table(as.matrix(mtx))
    mtx_dt[, "id"] <- as.character(add.id$V1)
    mtx_dt[, "ensembl.id"] <- as.character(add.id$V2)
    data.table::setkey(mtx_dt, ensembl.id, id)
    dt_ratio <- mtx_dt[, lapply(.SD, ratio), by = ensembl.id]
    DT.for.set.sqln(dt_ratio)
    RatioBysum <-
      Matrix::Matrix(as.matrix(dt_ratio[, -c("id", "ensembl.id")]), sparse = TRUE)
    rownames(RatioBysum) <- dt_ratio$id
    return(RatioBysum)
  }
  #split data,memory  would be 20% for each data when convert to data.table if set ncells =40000
  PAratio.mtx <- function(mtx,
                          add.id,
                          ncells = 40000,
                          cores = 1) {
    if (identical(as.character(add.id$V1), rownames(mtx))) {
      add.id <- add.id
    } else{
      rownames(add.id) <- add.id$V1
      add.id <- add.id[rownames(mtx),]
    }
    if (ncol(mtx) < ncells) {
      PAratio.mtx <- PAratio(mtx)
    } else {
      chunkwise <- seq(1, ncol(mtx), ncells)
      chunkwise.mtx <- list()
      for (i in 1:length(chunkwise)) {
        n <- chunkwise[i]
        N <- n + ncells - 1
        if (n < max(chunkwise)) {
          chunked.data <- mtx[, n:N]
          chunkwise.mtx[[i]] <- chunked.data
        } else{
          chunkwise.mtx[[i]] <- mtx[, n:ncol(mtx)]
        }
      }
      PAratio.mtx <-
        do.call(cbind,
                parallel::mclapply(chunkwise.mtx,
                                   function(mtx) {
                                     mtx_dt <- data.table::as.data.table(as.matrix(mtx))
                                     mtx_dt[, "id"] <-
                                       as.character(add.id$V1)
                                     mtx_dt[, "ensembl.id"] <-
                                       as.character(add.id$V2)
                                     data.table::setkey(mtx_dt, ensembl.id, id)
                                     dt_ratio <-
                                       mtx_dt[, lapply(.SD, ratio), by = ensembl.id]
                                     DT.for.set.sqln(dt_ratio)
                                     RatioBysum <-
                                       Matrix::Matrix(as.matrix(dt_ratio[, -c("id", "ensembl.id")]), sparse = TRUE)
                                     rownames(RatioBysum) <-
                                       dt_ratio$id
                                     return(RatioBysum)
                                   },
                                   mc.cores = cores))
    }
    return(PAratio.mtx)
  }

  Ratiobysum.matrix <-
    PAratio.mtx(pa.geneid.matrix,
                annot,
                ncells = chunk,
                cores = cores)
  obj[[assay]] <- CreateAssayObject(Ratiobysum.matrix)
  obj
}

#'@title Psi calculation.
#'@description Psi calculation.
#'@param obj A Seurat object include `apa` assay.
#'@param annot Annotation dataframe.
#'@param chunk Cell number per chunk for reducing memory usage.
#'@param cores The num of cpu for calculating psi.
#'@example
#'@example
#'@example
#'\dontrun{annot_info <-
#'   annotate_from_gtf('/home/zhouran/data/tmp/annotate_test/chr1.gtf.gz',
#'                     'Mm10',
#'                    cores = 10)
#'}
#'
#'@import data.table
#'@export

psiCate <- function(obj,
                    annot,
                    assay = 'apapsi') {
  # classification function
  ratio_classify <- function(psi.table) {
    fraction.psi <- apply(psi.table, 1, function(x) {
      ### calculate the observed statistical values
      out <- c(
        table(is.na(x))["TRUE"] / length(x),
        ### number of NAs
        table(x == 0)["TRUE"] / length(x),
        ### psi = 0
        table(x == 1)["TRUE"] / length(x),
        ### psi = 1
        table(x > 0 & x < 1)["TRUE"] / length(x),
        length(na.omit(x)),
        ### number of nonNA
        mean(x, na.rm = T),
        median(x, na.rm = T),
        var(x, na.rm = T),
        ### variation
        sd(x, na.rm = T),
        ### SD
        var(x, na.rm = T) / mean(x, na.rm = T) ### CV
      )
      out[is.na(out)] <- 0
      return(round(out, 3))
    })
    t(fraction.psi) -> a
    colnames(a) <- c(
      "frac.na",
      "frac.0",
      "frac.1",
      "frac.0-1",
      'num.nonNA',
      "mean_psi",
      "median",
      'Var_psi',
      'SD_psi',
      'CV_psi'
    )
    return(as.data.frame(a))
  }
  ratio.fraction <- function(psi.table) {
    fraction.psi <- apply(psi.table, 1, function(x) {
      #print(x)
      ### calculate the observed statistical values
      out <- c(
        table(is.na(x))["TRUE"] / length(x),
        ### number of NAs
        table(x == 0)["TRUE"] / length(x),
        ### psi = 0
        table(x == 1)["TRUE"] / length(x),
        ### psi = 1
        table(x > 0 & x < 1)["TRUE"] / length(x),
        length(na.omit(x)),
        ### number of nonNA
        mean(x, na.rm = T),
        median(x, na.rm = T),
        var(x, na.rm = T),
        ### variation
        sd(x, na.rm = T),
        ### SD
        var(x, na.rm = T) / mean(x, na.rm = T) ### CV
      )
      out[is.na(out)] <- 0
      return(round(out, 3))
    })
    t(fraction.psi) -> a
    colnames(a) <- c(
      "frac.na",
      "frac.0",
      "frac.1",
      "frac.0-1",
      'num.nonNA',
      "mean_psi",
      "median",
      'Var_psi',
      'SD_psi',
      'CV_psi'
    )
    return(as.data.frame(a))
  }

  psi_mat <- Seurat::GetAssayData(obj, assay = assay)
  add.id <- annot[rownames(psi_mat), ]

  tissue.ratio <-
    as.data.frame(as.matrix(psi_mat)) %>%
    dplyr::mutate(ensembl.id = as.character(add.id$annot.gene_id),
                  id = as.character(add.id$pa_site)) %>%
    data.table::as.data.table() %>%
    data.table::setkey(ensembl.id, id)


  tissue.data <-
    tissue.ratio[, lapply(.SD, function(x) {
      if (sum(x, na.rm = T) == 0) {
        x[x == 0] <- NA
      } else{
        x = x
      }
      x
    }), by = 'ensembl.id', .SDcols = !"id"]

  tissue.data[["id"]] <- tissue.ratio[["id"]]

  psi_mat <- psi_mat[Matrix::rowSums(psi_mat) > 0, ]

  tissue.use <-
    tissue.data[tissue.data$id %in% rownames(psi_mat), ]

  tissue.fractions <-
    ratio.fraction(tissue.use[, -c("id", "ensembl.id")])
  tissue.fractions$pa = tissue.use$id
  tissue.fractions <- data.table::as.data.table(tissue.fractions)

  arcsin_trans = t(apply(tissue.use[, -c("id", "ensembl.id")], 1, function(x) {
    asin(sqrt(x))
  }))

  tissue.fractions$Var_expect <-
    apply(arcsin_trans, 1, var, na.rm = T)
  tissue.fractions$SD_expect <-
    apply(arcsin_trans, 1, sd, na.rm = T)
  tissue.fractions$mean_expect <-
    apply(arcsin_trans, 1, mean, na.rm = T)

  tissue.fractions <-
    transform(tissue.fractions, Var_deviation = Var_expect - Var_psi)

  ### using the interquartile range (IQR) to define the category

  quantile(tissue.fractions$Var_deviation, na.rm = T)
  deviation_1st <-
    quantile(tissue.fractions$Var_deviation, na.rm = T)[3]
  deviation_3rd <-
    quantile(tissue.fractions$Var_deviation, na.rm = T)[4]


  tissue.fractions$PA_category <- 'Multimodal'
  tissue.fractions[mean_psi == 0]$PA_category <- 'NonExpr'
  tissue.fractions[mean_psi == 1]$PA_category <- 'NonAPA'
  tissue.fractions[mean_psi <= 0.2 &
                     mean_psi > 0]$PA_category <- 'L_shape'
  tissue.fractions[mean_psi >= 0.8 &
                     mean_psi < 1]$PA_category <- 'J_shape'
  tissue.fractions[mean_psi > 0.2 &
                     mean_psi < 0.8 &
                     Var_deviation >= deviation_3rd]$PA_category <-
    'OverDispersed'
  tissue.fractions[mean_psi > 0.2 &
                     mean_psi < 0.8 &
                     Var_deviation <= deviation_1st]$PA_category <-
    'UnderDispersed'

  tissue.fractions$PA_category <-
    as.character(tissue.fractions$PA_category)

  return(tissue.fractions)
}

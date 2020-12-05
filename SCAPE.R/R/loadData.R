library(data.table)
library(Seurat)

#'@title Loading datasets
#'@description loading apa matrix into a Seurat object.
#'@param fileList A named vector of apa matrix file.
#'@param collapsePa The output of `group_pa.py`.
#'@param matrix Return matrix or Seurat object.
#'@param cores The num of cpu for processing data.
#'@example
#'\dontrun{Not Run
#' fileList <- c(
#'  "Brain_1"="Brain1/pasite.csv.gz",
#'  "Brain_2"="Brain2/pasite.csv.gz",
#'  "BoneMarrow_1"="BoneMarrow1/pasite.csv.gz",
#'  "BoneMarrow_2"="BoneMarrow2/pasite.csv.gz",
#'  "BoneMarrow_3"="BoneMarrow3/pasite.csv.gz"
#'  )
#' dat <- loadData(fileList, 'collapse_ps.tsv.gz', cores=4)
#'}
#'
#'@export
loadData <- function(fileList,
                     collapsePa,
                     matrix = FALSE,
                     cores = 1) {
  pa_group_info <- data.table::fread(collapsePa)
  pa_group_info$pa <-
    paste(pa_group_info$chrom,
          pa_group_info$collapse_pa,
          pa_group_info$strand,
          sep = ':')

  objs <- parallel::mclapply(names(fileList), function(fileLabel) {
    message(fileLabel)


    # Load pA matrix information
    tmp <- data.table::fread(fileList[[fileLabel]])
    mtx <- as.matrix(tmp[,-1])
    rownames(mtx) <- tmp[['V1']]
    rm(tmp)
    # Mapping pA ID into collapse ID
    pa <- data.frame(V1 = rownames(mtx), stringsAsFactors = F)
    pa$V2 <-
      plyr::mapvalues(
        from = pa_group_info$pa_site,
        to = pa_group_info$pa,
        x = pa$V1,
        warn_missing = F
      )
    pa_dup <- pa[pa$V2 %in% pa$V2[duplicated(pa$V2)],]

    # Create new pA matrix
    mtx_dup <-
      do.call(rbind, lapply(split(pa_dup$V1, pa_dup$V2), function(x) {
        colSums(mtx[x, ])
      }))

    mtx <- rbind(mtx[setdiff(pa$V1, pa_dup$V1), ], mtx_dup)

    rownames(mtx) <-
      plyr::mapvalues(
        from = pa$V1,
        to = pa$V2,
        x = rownames(mtx),
        warn_missing = F
      )
    colnames(mtx) <- paste(fileLabel, colnames(mtx), sep = '.')

    dat <-
      Seurat::CreateSeuratObject(mtx, project = fileLabel, names.delim = '[.]')
    return(dat)
  }, mc.cores = cores)

  objs <- Reduce(merge, objs)
  if (isTRUE(matrix)) {
    return(Seurat::GetAssayData(objs, 'counts'))
  }
  return(objs)
}

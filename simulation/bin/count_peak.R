CountPeaks_fake <- function(peak.sites.file, 
                       gtf.file, 
                       bamfile, 
                       whitelist.file, 
                       output.dir, 
                       countUMI=TRUE,
			                 ncores = 1, 
			                 chr.names = NULL, 
			                 filter.chr = FALSE, 
			                 CBtag='CB', 
			                 UMItag='UB') 
{

  lock <- tempfile()
  whitelist.bc <- read.table(whitelist.file, stringsAsFactors = FALSE)
  whitelist.bc <- whitelist.bc[,1]
  n.bcs <- length(whitelist.bc)
  message("There are ", n.bcs, " whitelist barcodes.")

  n.columns <- n.bcs + 1

  # read in gene reference
  genes.ref <- make_reference(gtf_file = gtf.file, chr.names = chr.names, filter.chr = filter.chr)
  chr.names <- as.character(unique(genes.ref$chr))
  n.genes <- nrow(genes.ref)

  peak.sites <- read.table(peak.sites.file, header = T, sep = "\t",
                            stringsAsFactors = FALSE)

  # Count the peaks
  n.total.sites <- nrow(peak.sites)
  message("There are ", n.total.sites, "  sites")
  message("Doing counting for each site...")

  # Set up multiple workers
  system.name <- Sys.info()['sysname']
  new_cl <- FALSE
  if (system.name == "Windows") {
    new_cl <- TRUE
    cluster <- parallel::makePSOCKcluster(rep("localhost", ncores))
    doParallel::registerDoParallel(cluster)
  } else {
    doParallel::registerDoParallel(cores=ncores)
  }
  #print(chr.names)
  mat.to.write <- foreach::foreach(each.chr = chr.names, .combine = 'rbind', .packages=c("magrittr")) %dopar% {
#  mat.to.write <- foreach::foreach(each.chr = chr.names, .combine = 'rbind', .packages=c("magrittr")) %do% {
      mat.per.chr <- c()
      message("Processing chr: ", each.chr)
      

      for(strand in c(1, -1) ) {
      message(" and strand ", strand)
      isMinusStrand <- if(strand==1) FALSE else TRUE
     
      peak.sites.chr <- dplyr::filter(peak.sites, Chr == each.chr & Strand == strand) %>%
                           dplyr::select(Gene, Chr, Fit.start, Fit.end, Strand)

      peak.sites.chr$Fit.start <- as.integer(peak.sites.chr$Fit.start)
      peak.sites.chr$Fit.end <- as.integer(peak.sites.chr$Fit.end)
      peak.sites.chr <- dplyr::filter(peak.sites.chr, Fit.start < Fit.end)

      # If there are no sites in this range, then just keep going
      if(nrow(peak.sites.chr) == 0) {
	      next
      }

      isMinusStrand <- if(strand==1) FALSE else TRUE
      which <- GenomicRanges::GRanges(seqnames = each.chr, ranges = IRanges::IRanges(1, max(peak.sites.chr$Fit.end) ))

#      param <- Rsamtools::ScanBamParam(tag=c("CB", "UB"),     # CBtag='CB', UMItag='UB')
      param <- Rsamtools::ScanBamParam(tag=c(CBtag, UMItag),
                            which = which,
                            flag=Rsamtools::scanBamFlag(isMinusStrand=isMinusStrand))

      aln <- GenomicAlignments::readGAlignments(bamfile, param=param)

#      nobarcodes <- which(is.na(GenomicRanges::mcols(aln)$CB))
#      noUMI <- which(is.na(GenomicRanges::mcols(aln)$UB))
      GenomicRanges::mcols(aln)[CBtag] <- "AAA"
      nobarcodes <- which(unlist(is.na(GenomicRanges::mcols(aln)[CBtag])))
      GenomicRanges::mcols(aln)[UMItag] <- "BB"
      noUMI <- which(unlist(is.na(GenomicRanges::mcols(aln)[UMItag])))
      
      
      to.remove <- dplyr::union(nobarcodes, noUMI)
      if (length(to.remove) > 0) {
        aln <- aln[-to.remove]
      }
#      whitelist.pos <- which(GenomicRanges::mcols(aln)$CB %in% whitelist.bc)
      whitelist.pos <- which(unlist(GenomicRanges::mcols(aln)[CBtag]) %in% whitelist.bc)
      aln <- aln[whitelist.pos]

      # For de-duplicating UMIs, let's just remove a random read
      # when there is a duplicate
      if(countUMI) {
 #        GenomicRanges::mcols(aln)$CB_UB <- paste0(GenomicRanges::mcols(aln)$CB, "_", GenomicRanges::mcols(aln)$UB)
        GenomicRanges::mcols(aln)$CB_UB <- paste0(unlist(GenomicRanges::mcols(aln)[CBtag]), 
                                                  "_", unlist(GenomicRanges::mcols(aln)[UMItag]))
         uniqUMIs <- which(!duplicated(GenomicRanges::mcols(aln)$CB_UB))
         aln <- aln[uniqUMIs]
      }

#      aln <- GenomicRanges::split(aln, GenomicRanges::mcols(aln)$CB)
      aln <- GenomicRanges::split(aln, unlist(GenomicRanges::mcols(aln)[CBtag]))

      polyA.GR <- GenomicRanges::GRanges(seqnames = peak.sites.chr$Chr,
                          IRanges::IRanges(start = peak.sites.chr$Fit.start,
                                  end = as.integer(peak.sites.chr$Fit.end)))
      n.polyA <- length(polyA.GR)
      barcodes.gene <- names(aln)
      res <- sapply(barcodes.gene, function(x) GenomicRanges::countOverlaps(polyA.GR, aln[[x]]))

      # Reorder the columns of the res matrix to match the whitelist barcodes
      res.mat <- matrix(0L, nrow = n.polyA, ncol = n.bcs)
      res.mat[,match(barcodes.gene, whitelist.bc)] <- res

      # Return a sparse matrix
      mat.per.strand <- Matrix::Matrix(res.mat, sparse = TRUE)
      #mat.to.write <- matrix(0L, nrow = n.polyA, ncol = n.bcs)
      #mat.to.write[,match(barcodes.gene, whitelist.bc)] <- res
      polyA.ids <- paste0(peak.sites.chr$Gene, ":", peak.sites.chr$Chr, ":", peak.sites.chr$Fit.start,
                          "-", peak.sites.chr$Fit.end, ":", peak.sites.chr$Strand )
      rownames(mat.per.strand) <- polyA.ids


      # Need to combine the two matrices from each strand
      if(is.null(mat.per.chr)) {
         mat.per.chr <- mat.per.strand
      } else {
         mat.per.chr <- rbind(mat.per.chr, mat.per.strand)
      }
    } # Loop for strand

    # Return sparse matrix for each chromosome for combining across all threads
    return(mat.per.chr)
  } # Loop for chr
  
  if (new_cl) { ## Shut down cluster if on Windows
    ## stop cluster
    parallel::stopCluster(cluster)
  }

  if (!dir.exists(output.dir)){
    dir.create(output.dir)
  }
  Matrix::writeMM(mat.to.write, file = paste0(output.dir, "/matrix.mtx"))
  write.table(whitelist.bc, file = paste0(output.dir, "/barcodes.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(rownames(mat.to.write), file = paste0(output.dir, "/sitenames.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  ## Compress the output files
  R.utils::gzip(paste0(output.dir, "/matrix.mtx"), overwrite = TRUE)
  R.utils::gzip(paste0(output.dir, "/barcodes.tsv"), overwrite = TRUE)
  R.utils::gzip(paste0(output.dir, "/sitenames.tsv"), overwrite = TRUE)

} # End function

environment(CountPeaks_fake) <- asNamespace('Sierra')
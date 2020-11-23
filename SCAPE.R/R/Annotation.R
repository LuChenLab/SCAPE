# library(GenomicRanges)
# library(GenomicFeatures)
# library(annotatr)
# library(GenomicAlignments)


annotate_from_gtf <-
  function (gtfFile,
            genome_ver = c('Mm10', 'Hg38'),
            cores = 1) {
    genome_ver <- match.arg(genome_ver)
    species_map <- c(
      'Hg38' = 'Hg',
      'Hg19' = 'Hg',
      'Mm10' = 'Mm',
      'Mm9' = 'Mm'
    )

    if (species_map[[genome_ver]] == 'Mm') {
      require(org.Mm.eg.db)
      x_map = get(sprintf('org.Mm.egSYMBOL', genome_ver))
    } else if (species_map[[genome_ver]] == 'Hg') {
      require(org.Hs.eg.db)
      x_map = get(sprintf('org.Hs.egSYMBOL', genome_ver))
    } else {
      stop("Can't recognize the genome version!")

    }

    txdb <- GenomicFeatures::makeTxDbFromGFF(gtfFile, format = "gtf")
    ebg <- GenomicFeatures::transcriptsBy(txdb, by = "gene")


    mapped_genes = AnnotationDbi::mappedkeys(x_map)
    eg2symbol = as.data.frame(x_map[mapped_genes])
    eg2symbol$ensemble <-
      AnnotationDbi::mapIds(
        org.Mm.eg.db,
        keys = eg2symbol$gene_id,
        keytype = "ENTREZID",
        column = "ENSEMBL"
      )

    tx_gr = GenomicFeatures::transcripts(txdb, columns = c('TXID', 'GENEID', 'TXNAME'))

    id_maps = AnnotationDbi::select(
      txdb,
      keys = as.character(GenomicRanges::mcols(tx_gr)$TXID),
      columns = c('TXNAME', 'GENEID'),
      keytype = 'TXID'
    )

    # GenomicFeatures::exonsBy(txdb, by = 'tx', use.names = TRUE)
    exons_grl = GenomicFeatures::exonsBy(txdb, by = 'tx', use.names = TRUE)
    # Create Rle of the tx_names
    exons_txname_rle = S4Vectors::Rle(names(exons_grl), S4Vectors::elementNROWS(exons_grl))
    exons_txname_vec = as.character(exons_txname_rle)
    # Unlist and add the tx_names
    exons_gr = unlist(exons_grl, use.names = FALSE)
    GenomicRanges::mcols(exons_gr)$tx_name = exons_txname_vec
    # Add Entrez ID, symbol, and type
    GenomicRanges::mcols(exons_gr)$gene_id = id_maps[match(GenomicRanges::mcols(exons_gr)$tx_name, id_maps$TXNAME),
                                                     'GENEID']
    GenomicRanges::mcols(exons_gr)$symbol = eg2symbol[match(GenomicRanges::mcols(exons_gr)$gene_id,
                                                            eg2symbol$ensemble), 'symbol']
    GenomicRanges::mcols(exons_gr)$entrez = eg2symbol[match(GenomicRanges::mcols(exons_gr)$gene_id,
                                                            eg2symbol$ensemble), 'gene_id']
    # GenomicRanges::mcols(exons_gr)$type = sprintf('%s_genes_exons', )
    GenomicRanges::mcols(exons_gr)$id = paste0('Exon:ExonRank', exons_gr$exon_rank)


    x <- lapply(exons_grl, function(x) {
      exonrank <- length(x$exon_rank)
      if (as.character(strand(x)) == "+") {
        x[exonrank]
      } else{
        x[1]
      }
    })

    lastexons_grl <- GRangesList(x)

    lastexon_txname_rle = S4Vectors::Rle(names(lastexons_grl),
                                         S4Vectors::elementNROWS(lastexons_grl))
    lastexons_txname_vec = as.character(lastexon_txname_rle)

    lastexons_gr = unlist(lastexons_grl, use.names = FALSE)
    GenomicRanges::mcols(lastexons_gr)$tx_name = lastexons_txname_vec
    # Add Entrez ID, symbol, and type
    GenomicRanges::mcols(lastexons_gr)$gene_id = id_maps[match(GenomicRanges::mcols(lastexons_gr)$tx_name,
                                                               id_maps$TXNAME),
                                                         'GENEID']
    GenomicRanges::mcols(lastexons_gr)$symbol = eg2symbol[match(GenomicRanges::mcols(lastexons_gr)$gene_id,
                                                                eg2symbol$ensemble), 'symbol']
    GenomicRanges::mcols(lastexons_gr)$entrez = eg2symbol[match(GenomicRanges::mcols(lastexons_gr)$gene_id,
                                                                eg2symbol$ensemble), 'gene_id']
    # GenomicRanges::mcols(exons_gr)$type = sprintf('%s_genes_exons', )
    GenomicRanges::mcols(lastexons_gr)$id = paste0('lastexon:ExonRank', lastexons_gr$exon_rank)

    pos = lastexons_gr[strand(lastexons_gr) == '+', ]
    neg = lastexons_gr[strand(lastexons_gr) == '-', ]
    start(pos) <- end(pos) + 1
    end(pos) <- end(pos) + 1000
    end(neg) <- start(neg) - 1
    start(neg) <- start(neg) - 1000

    lastexons1k_gr_ <- c(pos, neg)
    GenomicRanges::mcols(lastexons1k_gr_)$id = paste0('LastExon1Kb:ExonRank', lastexons_gr$exon_rank)


    cds_grl = GenomicFeatures::cdsBy(txdb, by = 'tx', use.names = TRUE)
    cds_txname_rle = S4Vectors::Rle(names(cds_grl), S4Vectors::elementNROWS(cds_grl))
    cds_txname_vec = as.character(cds_txname_rle)
    # Unlist and add the tx_names
    cds_gr = unlist(cds_grl, use.names = FALSE)
    GenomicRanges::mcols(cds_gr)$tx_name = cds_txname_vec
    # Add Entrez ID, symbol, and type
    GenomicRanges::mcols(cds_gr)$gene_id = id_maps[match(GenomicRanges::mcols(cds_gr)$tx_name, id_maps$TXNAME), 'GENEID']
    GenomicRanges::mcols(cds_gr)$symbol = eg2symbol[match(GenomicRanges::mcols(cds_gr)$gene_id, eg2symbol$ensemble), 'symbol']
    GenomicRanges::mcols(cds_gr)$entrez = eg2symbol[match(GenomicRanges::mcols(cds_gr)$gene_id, eg2symbol$ensemble), 'gene_id']
    GenomicRanges::mcols(cds_gr)$type = sprintf('CDS')
    GenomicRanges::mcols(cds_gr)$id = paste0('CDS:ExonRank', cds_gr$exon_rank)


    fiveUTRs_grl = GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE)
    # Create Rle of the tx_names
    fiveUTRs_txname_rle = S4Vectors::Rle(names(fiveUTRs_grl),
                                         S4Vectors::elementNROWS(fiveUTRs_grl))
    fiveUTRs_txname_vec = as.character(fiveUTRs_txname_rle)
    # Unlist and add the tx_names
    fiveUTRs_gr = unlist(fiveUTRs_grl, use.names = FALSE)
    GenomicRanges::mcols(fiveUTRs_gr)$tx_name = fiveUTRs_txname_vec
    # Add Entrez ID, symbol, and type
    # NOTE: here we match on the tx_name because the tx_id is not given
    GenomicRanges::mcols(fiveUTRs_gr)$gene_id = id_maps[match(GenomicRanges::mcols(fiveUTRs_gr)$tx_name, id_maps$TXNAME), 'GENEID']
    GenomicRanges::mcols(fiveUTRs_gr)$symbol = eg2symbol[match(GenomicRanges::mcols(fiveUTRs_gr)$gene_id,
                                                               eg2symbol$ensemble), 'symbol']
    GenomicRanges::mcols(fiveUTRs_gr)$entrez = eg2symbol[match(GenomicRanges::mcols(fiveUTRs_gr)$gene_id,
                                                               eg2symbol$ensemble), 'gene_id']

    GenomicRanges::mcols(fiveUTRs_gr)$type = sprintf('5UTRs')
    GenomicRanges::mcols(fiveUTRs_gr)$id = paste0('5UTR:ExonRank', fiveUTRs_gr$exon_rank)


    threeUTRs_grl = GenomicFeatures::threeUTRsByTranscript(txdb, use.names = TRUE)
    # Create Rle of the tx_names
    threeUTRs_txname_rle = S4Vectors::Rle(names(threeUTRs_grl),
                                          S4Vectors::elementNROWS(threeUTRs_grl))
    threeUTRs_txname_vec = as.character(threeUTRs_txname_rle)
    # Unlist and add the tx_names
    threeUTRs_gr = unlist(threeUTRs_grl, use.names = FALSE)
    GenomicRanges::mcols(threeUTRs_gr)$tx_name = threeUTRs_txname_vec
    # NOTE: here we match on the tx_name because the tx_id is not given
    GenomicRanges::mcols(threeUTRs_gr)$gene_id = id_maps[match(GenomicRanges::mcols(threeUTRs_gr)$tx_name,
                                                               id_maps$TXNAME), 'GENEID']
    GenomicRanges::mcols(threeUTRs_gr)$symbol = eg2symbol[match(GenomicRanges::mcols(threeUTRs_gr)$gene_id,
                                                                eg2symbol$ensemble), 'symbol']
    GenomicRanges::mcols(threeUTRs_gr)$entrez = eg2symbol[match(GenomicRanges::mcols(threeUTRs_gr)$gene_id,
                                                                eg2symbol$ensemble), 'gene_id']

    GenomicRanges::mcols(threeUTRs_gr)$id = sprintf('3UTRs')

    pos = threeUTRs_gr[strand(threeUTRs_gr) == '+', ]
    neg = threeUTRs_gr[strand(threeUTRs_gr) == '-', ]

    start(pos) <- end(pos) + 1
    end(pos) <- end(pos) + 1000

    end(neg) <- start(neg) - 1
    start(neg) <- start(neg) - 1000

    threeUTRs1Kb <- c(neg, pos)
    threeUTRs1Kb$id <- '3UTRs_1kb'

    ###for 2kb
    pos = threeUTRs1Kb[strand(threeUTRs1Kb) == '+', ]
    neg = threeUTRs1Kb[strand(threeUTRs1Kb) == '-', ]

    start(pos) <- end(pos) + 1
    end(pos) <- end(pos) + 1000

    end(neg) <- start(neg) - 1
    start(neg) <- start(neg) - 1000

    threeUTRs2Kb <- c(pos, neg)
    threeUTRs2Kb$id <- '3UTRs_2kb'


    introns_grl = GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE)
    # Create Rle of the tx_names
    introns_txname_rle = S4Vectors::Rle(names(introns_grl), S4Vectors::elementNROWS(introns_grl))
    introns_txname_vec = as.character(introns_txname_rle)
    # Unlist and add the tx_names
    introns_gr = unlist(introns_grl, use.names = FALSE)
    GenomicRanges::mcols(introns_gr)$tx_name = introns_txname_vec
    # NOTE: here we match on the tx_name because the tx_id is not given
    GenomicRanges::mcols(introns_gr)$gene_id = id_maps[match(GenomicRanges::mcols(introns_gr)$tx_name, id_maps$TXNAME), 'GENEID']
    GenomicRanges::mcols(introns_gr)$symbol = eg2symbol[match(GenomicRanges::mcols(introns_gr)$gene_id,
                                                              eg2symbol$ensemble), 'symbol']
    GenomicRanges::mcols(introns_gr)$entrez = eg2symbol[match(GenomicRanges::mcols(introns_gr)$gene_id,
                                                              eg2symbol$ensemble), 'gene_id']



    introns_gr <-
      parallel::mclapply(split(introns_gr, introns_gr$tx_name), function(x) {
        if (unique(strand(x)) == "-") {
          x$id <- paste('Intron:Rank', rev(seq(1:length(x))), sep = '')
        } else{
          x$id <- paste('Intron:Rank', seq(1:length(x)), sep = '')
        }
        x
      }, mc.cores = cores)
    names(introns_gr) <- NULL
    introns_gr <- unlist(as(introns_gr, "GRangesList"))




    promoters_gr = GenomicFeatures::promoters(txdb, upstream = 1000, downstream = 0)
    # Add Entrez ID, symbol, and type
    GenomicRanges::mcols(promoters_gr)$gene_id = id_maps[match(GenomicRanges::mcols(promoters_gr)$tx_id, id_maps$TXID), 'GENEID']
    GenomicRanges::mcols(promoters_gr)$symbol = eg2symbol[match(GenomicRanges::mcols(promoters_gr)$gene_id,
                                                                eg2symbol$gene_id), 'symbol']
    GenomicRanges::mcols(promoters_gr)$type = sprintf('%s_genes_promoters', genome_ver)
    GenomicRanges::mcols(promoters_gr)$id = paste0('promoter:', seq_along(promoters_gr))

    GenomicRanges::mcols(promoters_gr) = GenomicRanges::mcols(promoters_gr)[, c('id', 'tx_name', 'gene_id', 'symbol', 'type')]
    colnames(GenomicRanges::mcols(promoters_gr)) = c('id', 'tx_id', 'gene_id', 'symbol', 'type')


    onetofive_gr = GenomicRanges::flank(promoters_gr,
                                        width = 4000,
                                        start = TRUE,
                                        both = FALSE)
    onetofive_gr = GenomicRanges::trim(onetofive_gr)
    # Add Entrez ID, symbol, and type (all but type are inherited from promoters_gr)
    GenomicRanges::mcols(onetofive_gr)$id = paste0('1to5kb:', seq_along(onetofive_gr))
    GenomicRanges::mcols(onetofive_gr)$type = sprintf('%s_genes_1to5kb', genome_ver)


    genic_gr = c(GenomicRanges::granges(tx_gr))
    GenomicRanges::strand(genic_gr) = '*'
    intergenic_gr =  (genic_gr)
    intergenic_gr = GenomicRanges::gaps(intergenic_gr)

    # A quirk in gaps gives the entire + and - strand of a chromosome, ignore those
    intergenic_gr = intergenic_gr[GenomicRanges::strand(intergenic_gr) == '*']

    GenomicRanges::mcols(intergenic_gr)$tx_name = 'NA'
    GenomicRanges::mcols(intergenic_gr)$gene_id = 'NA'
    GenomicRanges::mcols(intergenic_gr)$symbol = 'NA'
    GenomicRanges::mcols(intergenic_gr)$entrez = 'NA'
    GenomicRanges::mcols(intergenic_gr)$id = 'INTERGENIC'

    use_ind <- c("tx_name", "gene_id", "symbol", "entrez", "id")
    ##for cds
    cds_gr <- cds_gr[, use_ind]

    ##for 3utrs
    threeUTRs_gr <- threeUTRs_gr[, use_ind]

    ##for 3utrs-1kb
    threeUTRs1Kb <- threeUTRs1Kb[, use_ind]

    ##for 3utrs-2kb
    threeUTRs2Kb <- threeUTRs2Kb[, use_ind]

    ##for 5UTR
    fiveUTRs_gr <- fiveUTRs_gr[, use_ind]

    ##fo intergenic
    intergenic_gr <- intergenic_gr[, use_ind]

    ##for intron
    introns_gr <- introns_gr[, use_ind]
    ## for exon
    exons_gr <- exons_gr[, use_ind]

    ## 2019-6-14 modified to remove 5'utr extend

    annotation_db <- c(
      introns_gr,
      cds_gr,
      exons_gr,
      threeUTRs1Kb,
      threeUTRs_gr,
      fiveUTRs_gr,
      threeUTRs2Kb,
      intergenic_gr,
      lastexons1k_gr_
    )

    return(annotation_db)
  }


#'@title Annotate pA site based gtf file.
#'@description Annotate pA site based gtf file.
#'@param pAsite A list of pA sites, `1:810052:+`
#'@param gtfFile The standard gtf file.
#'@param genome_ver The version of genome annotation. Hg38, Hg19, Mm10, Mm9.
#'@param annotLevels The priority of annotaiton.
#'@param cores The num of cpu for DE test.
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

AnnotationSite <-
  function(pAsite,
           gtfFile,
           genome_ver,
           annotLevels = c(
             "3UTRs",
             "Exon",
             "Intron",
             "CDS",
             "LastExon1Kb",
             "3UTRs_1kb",
             "3UTRs_2kb",
             "5UTR",
             "INTERGENIC"
           ),
           cores = 1) {
    pa_info <-
      as.data.frame(do.call(rbind, (strsplit(pAsite, split = ':'))))
    colnames(pa_info) <- c('chr', 'end', 'strand')
    pa_info$end <- as.numeric(as.character(pa_info$end))
    pa_info$start <- pa_info$end - 1
    pa_info$score <- '.'
    pa_info <-
      GenomicRanges::makeGRangesFromDataFrame(pa_info, keep.extra.columns = TRUE)
    gtf_prefix <- strsplit(basename(gtfFile), split = '[.]')[[1]][1]
    annot_db_file <-
      file.path(dirname(gtfFile), paste0(gtf_prefix, '.Rds'))

    if (file.exists(annot_db_file)) {
      annot_db <- readRDS(annot_db_file)
    } else {
      annot_db <-
        annotate_from_gtf(gtfFile = gtfFile,
                          genome_ver = genome_ver,
                          cores = cores)
      saveRDS(annot_db, file = annot_db_file)
    }

    annot_res <- annotatr::annotate_regions(
      regions = pa_info,
      annotations = annot_db,
      ignore.strand = FALSE,
      quiet = FALSE
    )
    annot_res <-
      as.data.frame(annot_res, row.names = 1:length(annot_res))
    inds <-
      c(
        "seqnames",
        "end",
        "strand",
        "annot.start",
        "annot.end",
        "annot.tx_name",
        "annot.gene_id",
        "annot.symbol",
        "annot.entrez",
        "annot.id"
      )

    annot_res <- annot_res[, inds]
    annot_res <-
      cbind(annot_res, as.data.frame(do.call(
        rbind, strsplit(annot_res$annot.id, split = ':')
      )))
    colnames(annot_res)[c(11, 12)] <- c('Type', 'Rank')
    annot_res <- annot_res[,-10]
    annot_res$Type <- factor(annot_res$Type, levels = annotLevels)
    annot_res_lst <-
      split(annot_res,
            paste(annot_res$seqnames, annot_res$end, annot_res$strand, sep = ':'))
    annot_res_lst <- lapply(annot_res_lst, function(x) {
      ind <- order(x$Type)[1]
      return(x[ind,])
    })
    annot_res_lst <- do.call(rbind, annot_res_lst)
    annot_res_lst$pa_site <- rownames(annot_res_lst)
    colnames(annot_res)[2] <- 'pa_loc'
    return(annot_res_lst)
  }



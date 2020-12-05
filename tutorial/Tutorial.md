Introduction
============

This is a quick walkthrough of the downstream analysis of **SCAPE**.
The bone marrow, brain, fetal liver and spleen datasets from Microwell-seq have been processed by SCAPE.

Load the gene expression matrix into Seurat
-------------------------------------------

Load gene expression and umap information from prepared datasets.

``` r
library(Seurat)
library(SCAPE)
library(magrittr)
# Load gene expression
gene_obj <-
  readRDS(system.file('extdata', 'gene_mtx.Rds', package = 'SCAPE'))

# Load pre-calculated umap information
umap_coord <-
  readRDS(system.file('extdata', 'umap.Rds', package = 'SCAPE'))

# Load cell type
cell_ident <-
  readRDS(system.file('extdata', 'cellIdent.Rds', package = 'SCAPE'))

# Create Seurat object
gene_obj <-
  Seurat::CreateSeuratObject(counts = gene_obj, names.delim = '[.]')
```

``` r
gene_obj[["percent.mt"]] <-
  PercentageFeatureSet(gene_obj, pattern = "^mt-")

VlnPlot(
  gene_obj,
  features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'),
  group.by = 'orig.ident',
  ncol = 3
)
```

![](https://github.com/zhou-ran/SCAPE/blob/main/tutorial/Tutorial_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
gene_obj %<>% NormalizeData %<>% ScaleData

gene_obj[['umap']] <- umap_coord
gene_obj[['cellIdent']] <-
  plyr::mapvalues(
    from = rownames(cell_ident),
    to = as.character(cell_ident$cell_annot_collapse),
    x = colnames(gene_obj)
  )

Idents(gene_obj) <- 'cellIdent'

UMAPPlot(gene_obj, label = T) + NoLegend()
```

![](https://github.com/zhou-ran/SCAPE/blob/main/tutorial/Tutorial_files/figure-markdown_github/unnamed-chunk-3-1.png)

Load the APA expression matrix into Seurat.
-------------------------------------------

``` r
dir_loc <- system.file("extdata", "", package = "SCAPE")
files <- list.files(dir_loc, full.names = T, recursive = T)

# selecet the expression file
exp_file <- grep('pasite.csv.gz', files, value = T)
names(exp_file) <- basename(dirname(exp_file))

# load the collapse pA site file.
# generate from `script/group_pa.py`
collapse_pa <-
  system.file("extdata",
              "collapse_pa.tsv.gz",
              package = "SCAPE")

pa_mtx <- loadData(
  fileList = exp_file,
  collapsePa = collapse_pa,
  matrix = TRUE,
  cores = 8
)
```

Load pa matrix into Seurat object
---------------------------------

``` r
# Only these pA sites whcih expressed in more than 50 cell were kept.
binary_filter <- Matrix::rowSums(+pa_mtx)
pa_mtx <- pa_mtx[binary_filter > 50, ]

gene_obj[['apa']] <- CreateAssayObject(pa_mtx[, colnames(gene_obj)])
gene_obj <- NormalizeData(gene_obj, assay = 'apa')
gene_obj <- ScaleData(gene_obj, assay = 'apa')
```

``` r
VlnPlot(
  gene_obj,
  features = c('nFeature_RNA', 'nCount_RNA'),
  group.by = 'orig.ident',
  assay = 'apa',
  ncol = 2
)
```

![](https://github.com/zhou-ran/SCAPE/blob/main/tutorial/Tutorial_files/figure-markdown_github/unnamed-chunk-6-1.png)

Annotation of pA
----------------

``` r
gtf_file <-
  system.file('extdata', 'GRCm38.p5.genes.gtf.gz', package = 'SCAPE')

# It will consume a lot of time if it is the first time to annotate.

annot_info <-
  AnnotationSite(rownames(GetAssayData(gene_obj, assay = 'apa')),
                 gtf_file,
                 'Mm10',
                 cores = 10)
```

Find differential apa events
----------------------------

``` r
de_res <- SCAPE::FindDE(
  gene_obj,
  idents.1 = 'Erythroblast',
  idents.2 = 'Myelinating_oligodendrocyte',
  annot = annot_info,
  assay = "apa",
  cores = 20
)
```

Visualization of APA events

``` r
gene_to_check <- c('Bzw1')
pa <- c('1:58407357:+', '1:58405443:+')

FeaturePlot(gene_obj,
                 c(gene_to_check, pa),
                 order = T,
                 ncol = 3)
```

![](https://github.com/zhou-ran/SCAPE/blob/main/tutorial/Tutorial_files/figure-markdown_github/unnamed-chunk-9-1.png)

----------------------------
![](https://github.com/zhou-ran/SCAPE/blob/main/tutorial/Tutorial_files/figure-markdown_github/Bzw1.png)

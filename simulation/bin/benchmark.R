
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(MASS)
library(LSD)

generate_evaluation <- function(dat, truth) {
  gene_id <- names(truth)
  gene_id_use <- intersect(gene_id, dat$gene_id)
  
  precision <- lapply(seq(30, 150, 20), function(trial) {
    padistrance <- lapply(gene_id_use, function(x) {
      inferred <- dat[dat$gene_id == x, ]
      inferred_pa <- as.numeric(inferred$pa)
      ground_truth <- truth[[x]]
      pa_site <- lapply(inferred_pa, function(alpha_x) {
        if (any(abs(alpha_x - ground_truth$pasite) < trial)) {
          return(c(alpha_x, ground_truth$pasite[which.min(abs(alpha_x - ground_truth$pasite))]))
        } else {
          return(NULL)
        }
      })
      pa_site <-
        data.frame(do.call(rbind, pa_site[unlist(lapply(pa_site, function(x) {
          !is.null(x)
        }))]))
      
      if (isTRUE(pa_site)) {
        colnames(pa_site) <- c('inferred', 'truth')
      }
      
      pa_site
    })
    
    padistrance <- do.call(rbind, padistrance)
    
    tmp <-
      dim(padistrance)[1] / sum(sapply(lapply(truth, dim), '[[', 1))
    tmp
  })
  precision <-
    data.table::data.table(cbind(
      cutoff = seq(30, 150, 20),
      precision = unlist(precision)
    ))
  
  check_status <- all(is.na(dat$pa_counts))
  dat <- split(dat, dat$gene_id)
  if (isTRUE(check_status)) {
    dat <- lapply(dat, function(x) {
      x$ws <- x$pa_weights / sum(x$pa_weights)
      x
    })
    
  } else {
    dat <- lapply(dat, function(x) {
      x$ws <- x$pa_counts / sum(x$pa_counts)
      x
    })
    
  }
  
  
  ws <- lapply(gene_id_use, function(x) {
    lapply(seq(30, 150, 20), function(trial) {
      inferred <- dat[[x]]
      inferred_pa <- as.numeric(inferred$pa)
      ws <- inferred$ws
      
      ground_truth <- truth[[x]]
      ground_ws <-
        ground_truth$coverage / sum(ground_truth$coverage)
      ws_info <- lapply(1:length(inferred_pa), function(inds) {
        alpha_x <- inferred_pa[inds]
        if (any(abs(alpha_x - ground_truth$pasite) < trial)) {
          return(c(ws[inds], ground_ws[which.min(abs(alpha_x - ground_truth$pasite))]))
        } else {
          return(NULL)
        }
      })
      
      ws_info <-
        data.frame(do.call(rbind, ws_info[unlist(lapply(ws_info, function(x) {
          !is.null(x)
        }))]))
      
      if (sum(dim(ws_info)) != 0) {
        colnames(ws_info) <- c('inferred', 'truth')
        ws_info$cut_off <- trial
      }
      ws_info
    }) -> ws_infos
    do.call(rbind, ws_infos)
  })
  
  ws <- data.table::data.table(do.call(rbind, ws))
  
  pa_found <- sum(sapply(lapply(dat, dim), '[[', 1))
  
  res <- list(precision = precision,
              ws = ws,
              pa_found = pa_found)
  return(res)
  
}

generate_plot_mat <- function(bench, truth) {
  if (is.null(names(bench))) {
    stop('No names found in bench dataset')
  }
  # for precision
  lapply(names(bench), function(x) {
    pa_found <- data.table(table(bench[[x]]$ws$cut_off))
    tmp <- bench[[x]][['precision']]
    tmp$tool <- x
    tmp$recall <-
      tmp$precision * pa_found$N / sum(sapply(lapply(truth, dim), '[[', 1))
    tmp$f_score <-
      2 * tmp$precision * tmp$recall / (tmp$precision + tmp$recall)
    tmp
  }) %>% do.call(what = rbind, args = .) -> precision
  
  # for ws
  lapply(names(bench), function(x) {
    tmp <- bench[[x]][['ws']]
    tmp$tool <- x
    tmp
  }) %>% do.call(what = rbind, args = .) -> ws
  
  list(prf = precision,
       ws = ws)
}

# ws_plot <- function(mtx, label) {
#   num_sam <- dim(mtx)[1]
#   if (isTRUE(num_sam < 200)) {
#     p <-
#       ggplot(mtx, aes(x = inferred, y = truth)) +
#       theme_bw() +
#       geom_point(alpha = 0.5) +
#       xlab(glue::glue("{label}, N={num_sam}"))  + ylab('Ground truth') +
#       xlim(0, 1) + ylim(0, 1) +
#       coord_equal()
#     return(p)
#   }
#   dens <- kde2d(mtx$inferred, mtx$truth)
  
#   # create a new data frame of that 2d density grid
#   # (needs checking that I haven't stuffed up the order here of z?)
#   gr <- data.frame(with(dens, expand.grid(x, y)), as.vector(dens$z))
#   names(gr) <- c("xgr", "ygr", "zgr")
  
#   # Fit a model
#   mod <- loess(zgr ~ xgr * ygr, data = gr)
  
#   # Apply the model to the original data to estimate density at that point
#   mtx$pointdens <-
#     predict(mod, newdata = data.frame(xgr = mtx$inferred, ygr = mtx$truth))
  
  
#   # Draw plot
#   p <-
#     ggplot(mtx, aes(x = inferred, y = truth, color = pointdens)) +
#     theme_bw() +
#     scale_colour_gradientn(colours = colorpalette('heat', 5)) + xlim(0, 1) + ylim(0, 1)
#   p <-
#     p + geom_density2d(color = 'black') + ggrastr::geom_point_rast(alpha =
#                                                                      0.5)
#   p <-
#     p + xlab(glue::glue("{label}, N={num_sam}")) + ylab('Ground truth') +  geom_abline(slope = 1,
#                                                                                        intercept = 0,
#                                                                                        linetype = 3) +
#     stat_cor(
#       data = mtx[, c('inferred', 'truth')],
#       method = "pearson",
#       size = 3,
#       col = "red",
#       label.y.npc = "top",
#       label.x.npc = "left",
#       aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))
#     ) + coord_equal()# + ggtitle(glue::glue("{label} inferred"))
#   return(p)
  
# }
ws_plot <- function(mtx, label) {
  num_sam <- dim(mtx)[1]
  # if (isTRUE(num_sam < 200)) {
  #   p <-
  #     ggplot(mtx, aes(x = inferred, y = truth)) +
  #     theme_bw() +
  #     geom_point(alpha = 0.5) +
  #     xlab(glue::glue("{label}, N={num_sam}"))  + ylab('Ground truth') +
  #     xlim(0, 1) + ylim(0, 1) +
  #     coord_equal()
  #   return(p)
  # }
  h <- c(
    ifelse(
      bandwidth.nrd(mtx$inferred) == 0,
      0.1,
      bandwidth.nrd(mtx$inferred)
    ),
    ifelse(bandwidth.nrd(mtx$truth) == 0, 0.1, bandwidth.nrd(mtx$truth))
  )
  dens <- kde2d(mtx$inferred, mtx$truth,  h = h)
  
  # create a new data frame of that 2d density grid
  # (needs checking that I haven't stuffed up the order here of z?)
  gr <- data.frame(with(dens, expand.grid(x, y)), as.vector(dens$z))
  names(gr) <- c("xgr", "ygr", "zgr")
  
  # Fit a model
  mod <- loess(zgr ~ xgr * ygr, data = gr)
  
  # Apply the model to the original data to estimate density at that point
  mtx$pointdens <-
    predict(mod, newdata = data.frame(xgr = mtx$inferred, ygr = mtx$truth))
  
  
  # Draw plot
  p <-
    ggplot(mtx, aes(x = inferred, y = truth, color = pointdens)) +
    theme_bw() +
    scale_colour_gradientn(colours = colorpalette('heat', 5)) + xlim(0, 1) + ylim(0, 1)
  p <-
    p + geom_density2d(color = 'black', h = h) + ggrastr::geom_point_rast(alpha =
                                                                                     0.5)
  p <-
    p + xlab(glue::glue("{label}, N={num_sam}")) + ylab('Ground truth') +  geom_abline(slope = 1,
                                                                                       intercept = 0,
                                                                                       linetype = 3) +
    stat_cor(
      data = mtx[, c('inferred', 'truth')],
      method = "pearson",
      size = 3,
      col = "red",
      label.y.npc = "top",
      label.x.npc = "left",
      aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))
    ) + coord_equal()# + ggtitle(glue::glue("{label} inferred"))
  return(p)
  
}


list(
  MAAPER = 'MAAPER/MAAPER.Rds',
  scAPA = 'scAPA/scAPA.Rds',
  scAPAtrap = 'scAPAtrap/scAPAtrap.Rds',
  SCAPTURE = 'SCAPTURE/SCAPTURE.Rds',
  scDaPars = 'scDapars/scDapars.Rds',
  Sierra = 'Sierra/Sierra.Rds',
  polyApipe = 'polyApipe/polyApipe.Rds'
) -> dats


truth <- readRDS('ground_truth.Rds')
dats <- lapply(dats, readRDS)

bench <-
  lapply(names(dats), function(label) {
    message(label)
    dat <- dats[[label]]
    dd <- generate_evaluation(dat = dat, truth = truth)
  })
names(bench) <- names(dats)

bench_pro <- generate_plot_mat(bench, truth)

save(dats, bench, bench_pro, truth, file = 'benchmarks.RData')

dir.create('plot', showWarnings = F)
p1 <-
  ggplot(reshape2::melt(bench_pro$prf, id = c('tool', 'cutoff')),
         aes(x = cutoff, y = value, color = tool)) +
  geom_line(size = 1) +
  geom_point(size = 1) +
  theme_bw() +
  ylim(0, 1) +
  scale_x_continuous(breaks = seq(30, 150, 20)) +
  xlab('Cutoff(bp)') +
  Seurat:::FacetTheme() +
  facet_wrap(variable ~ ., ncol = 3)

lapply(split(bench_pro$ws, bench_pro$ws$tool), function(x){
  p <- ws_plot(mtx = x[x$cut_off == 50], label = unique(x$tool))
  p
}) -> p2

pdf.options(paper = 'a4')
pdf('plot/f_score.pdf', 10, 2.5)
p1
dev.off()

pdf('plot/ws_merge.pdf')
cowplot::plot_grid(plotlist = lapply(p2, function(x) {
  x + Seurat::NoLegend()
}),
nrow = 3,
ncol = 3)
dev.off()

pdf('plot/ws.pdf', 4, 4)
p2
dev.off()



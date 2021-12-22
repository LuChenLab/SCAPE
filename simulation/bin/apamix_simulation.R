# library(MCMCpack, lib.loc = '/mnt/data8/zhouran/lib/3.6')
############  auxiliary functions  ######################
# generate integers that are normally distributed, set to ymin or ymax if lower or higher
gen_norm_int = function(n, # generate n numbers
                        mu, # mean
                        sigma, # std
                        ymin=NA, # min
                        ymax=NA # max
){
  mu = round(mu)
  sigma = round(sigma)
  if(is.na(ymin)){
    ymin = mu - 3*sigma
  }
  if(is.na(ymax)){
    ymax = mu + 3*sigma
  }
  
  y = round(rnorm(n, mu, sigma))
  y[y<ymin] = ymin
  y[y>ymax] = ymax
  
  return(y)
}

# function for simulating APA reads
apa_sim = function(utr_len=1500,  # UTR length
                   
                   # parameters for polya length
                   max_polya_len=150, # maximum polya length
                   min_polya_len=5, # minimum polya length
                   empirical_polya_len_flag=TRUE,  # use empirical distribution or not
                   empirical_polya_len_arr=c(10, 30, 50, 70, 90, 110, 130),  # polya length arr
                   empirical_polya_len_pmf=c(309912, 4107929, 802856, 518229, 188316, 263208, 101), # polya length pmf
                   
                   # parameters for PA sites
                   mu_pa_sites=c(500, 1000), # mean of PA site locations
                   std_pa_sites=c(10, 30),
                   pa_ws=c(0.5, 0.5),    # weights for PA sites
                   noise_ws=0.1, # weights for noise component
                   
                   # parameters for fragment length
                   mu_frag = 300, # mean of fragment size
                   std_frag = 30, # standard deviation of fragment size
                   dispersion_rate = 1, # dispersion_rate x std_frag is actual std
                   
                   # parameters for read 1 length
                   mu_read_utr = 120, # mean of r1 length on UTR
                   std_read_utr = 10, # std of r1 length on UTR
                   
                   # parameters for read 2 length, sampled from a mixture of Gaussian and uniform distribution
                   mu_read_polya = 30, # mean of r2 length on polya
                   std_read_polya = 5, # std of r2 length of polya
                   unif_ws_read_polya = 0.1, # probability sampled from uniform distribution
                   
                   # parameter for captured polya length
                   map_rate = 0.1, # probability of r2 that can be mapped to UTR
                   
                   n_frag=2000 # number of fragments to be generated
                   
){
  
  ############  part 1:  PA sites  ######################
  n_pa_site = length(mu_pa_sites)
  
  ws = c(pa_ws, noise_ws)
  ws = ws/sum(ws)
  pa_ws = ws[1:n_pa_site]
  noise_ws = ws[n_pa_site + 1]
  
  ############  part 2:  fragments and reads ######################
  
  reads_boarder = round(n_frag * ws)
  reads_boarder = cumsum(reads_boarder)
  reads_boarder[length(reads_boarder)] = n_frag
  reads_boarder = c(0, reads_boarder)
  
  frag_label = rep(0, n_frag)
  
  # generate polya length
  if (empirical_polya_len_flag) {
    # 根据经验分布来随机生成polyA的长度
    # generate using empirical distribution, FLAM paper
    polya_len_res = ksmooth(
      x = empirical_polya_len_arr,
      y = empirical_polya_len_pmf,
      kernel = 'normal',
      bandwidth = 20,
      x.points = seq(min_polya_len, max_polya_len)
    )
    polya_part_len_arr = sample(polya_len_res$x,
                                n_frag,
                                replace = T,
                                prob = polya_len_res$y)  # s
  } else{
    # sample uniformly
    # 如果不同经验分布用uniform distribution
    polya_part_len_arr = sample(seq(min_polya_len, max_polya_len), n_frag, replace = T)
  }
  
  # generate fragment size
  frag_len_arr = gen_norm_int(n_frag, mu_frag, std_frag * dispersion_rate)  # theta - x + 1 + s
  
  # generate fragment start position and PA sites
  # 生成 x, start site
  read_utr_st_arr = rep(0, n_frag)    # x
  component_cnt_arr = rep(0, n_pa_site + 1) # counts of each component
  theta_arr = rep(NA, n_frag)   # PA site for reads
  for(i in seq(n_pa_site+1)){
    tmpinds = seq((reads_boarder[i] + 1), reads_boarder[i + 1])
    component_cnt_arr[i] = length(tmpinds)
    frag_label[tmpinds] = i
    if(i<=n_pa_site){
      theta_arr[tmpinds] = gen_norm_int(length(tmpinds), mu_pa_sites[i], std_pa_sites[i], 1, utr_len)
      read_utr_st_arr[tmpinds] = -frag_len_arr[tmpinds] + theta_arr[tmpinds] + 1 + polya_part_len_arr[tmpinds]
      if(any(read_utr_st_arr[tmpinds]<=0)){
        stop(paste0('Start of some fragments for pa site i=',i,' (mu=', mu_pa_sites[i], ' std=', std_pa_sites[i], ') is less than 0. Consider shorten fragment length or move pa site towards polya.\n'))
      }
    }else{
      read_utr_st_arr[tmpinds] = round( runif(length(tmpinds))*(utr_len-mu_read_utr-std_read_utr*3) )
    }
  }
  
  # generate read lenghts on UTR and polyA, trim read if longer than fragment
  # 生成 l
  utr_read_len_arr = gen_norm_int(n_frag, mu_read_utr, std_read_utr)   # l
  
  # r from Gaussian
  polya_read_len_arr = rep(0, n_frag)
  gaussian_flag = runif(n_frag) >= unif_ws_read_polya
  polya_read_len_arr[gaussian_flag] = gen_norm_int(sum(gaussian_flag), mu_read_polya, std_read_polya)   # r
  # r from uniform
  unif_inds = which(!gaussian_flag)
  for(i in unif_inds){
    polya_read_len_arr[i] = sample(polya_part_len_arr[i],1)
  }
  
  # 防止r2比对长度大于fragment length
  pa_inds = which(frag_label != (n_pa_site + 1))
  if(any(utr_read_len_arr[pa_inds]>frag_len_arr[pa_inds])){
    tmpinds0 = utr_read_len_arr[pa_inds]>frag_len_arr[pa_inds]
    tmpinds1 = pa_inds[tmpinds0]
    utr_read_len_arr[tmpinds1] = frag_len_arr[tmpinds1]
    warning(paste0('Length of utr reads longer than fragment length. Trimed to fit fragment length.\n',
                   'index=',paste(tmpinds1,collapse=" "),"\n"
    ))
  }
  
  # 防止polya length 长度大于fragment length
  if (any(polya_read_len_arr[pa_inds] > frag_len_arr[pa_inds])) {
    tmpinds0 = polya_read_len_arr[pa_inds] > frag_len_arr[pa_inds]
    tmpinds1 = pa_inds[tmpinds0]
    polya_read_len_arr[tmpinds1] = frag_len_arr[tmpinds1]
    warning(
      paste0(
        'Length of polya reads longer than fragment length. Trimed to fit fragment length.\n',
        'index=',
        paste(tmpinds1, collapse = " "),
        "\n"
      )
    )
  }
  
  noise_inds = which(frag_label == (n_pa_site + 1))
  s = rep(NA,n_frag)
  t = rep(NA,n_frag)
  
  # check if reads touch PA sites from utr side (APA component and noise component)
  # 2021-11-8 remove this part to generate R2 with polyA

  # tmpinds = read_utr_st_arr[pa_inds] + utr_read_len_arr[pa_inds] - 1 >= theta_arr[pa_inds]
  # if(isTRUE(any(tmpinds))){
  #   utr_read_len_arr[pa_inds[tmpinds]] = -read_utr_st_arr[pa_inds[tmpinds]]+1+theta_arr[pa_inds[tmpinds]]
  #   t[pa_inds[tmpinds]] = theta_arr[pa_inds[tmpinds]]
  # }
  # tmpinds = read_utr_st_arr[noise_inds] + utr_read_len_arr[noise_inds] - 1 >= utr_len
  # if(isTRUE(any(tmpinds))){
  #   utr_read_len_arr[noise_inds[tmpinds]] = -read_utr_st_arr[noise_inds[tmpinds]]+1+utr_len
  #   t[noise_inds[tmpinds]] = theta_arr[noise_inds[tmpinds]]
  # }
  
  # perform bernouli experiment to simulate the mapping process on the polya side
  r = pmin(polya_read_len_arr, polya_part_len_arr)
  
  bern_arr = runif(n_frag)
  map_flag_arr = bern_arr < map_rate
  map_inds = which(map_flag_arr)
  if(length(map_inds)>0){
    tmpinds = utr_read_len_arr[map_inds] >= polya_part_len_arr[map_inds]
    t[map_inds[tmpinds]]=theta_arr[map_inds[tmpinds]]
    s[map_inds[tmpinds]]=polya_part_len_arr[map_inds[tmpinds]]
    r[map_inds[tmpinds]] = polya_part_len_arr[map_inds[tmpinds]]
  }
  
  # generate inputs for apamix 
  L = utr_len
  LA = max_polya_len
  x = read_utr_st_arr
  l = utr_read_len_arr
  
  z = frag_label
  ws = component_cnt_arr
  ws = ws/sum(ws)
  
  return(list(L=L, LA=LA, x=x, l=l, r=r, s=s, t=t, 
              theta_arr=theta_arr, utr_read_len_arr=utr_read_len_arr,polya_part_len_arr=polya_part_len_arr,
              z=z, ws=ws,mu_pa_sites=mu_pa_sites,utr_read_len_arr=utr_read_len_arr,
              polya_read_len_arr=polya_read_len_arr,component_cnt_arr=component_cnt_arr,
              # input parameters
              utr_len=utr_len,  # UTR length
              
              # parameters for polya length
              max_polya_len=max_polya_len, # maximum polya length
              min_polya_len=min_polya_len, # minimum polya length
              empirical_polya_len_flag=empirical_polya_len_flag,  # use empirical distribution or not
              empirical_polya_len_arr=empirical_polya_len_arr,  # polya length arr
              empirical_polya_len_pmf=empirical_polya_len_pmf, # polya length pmf
              
              # parameters for PA sites
              mu_pa_sites=mu_pa_sites, # mean of PA site locations
              std_pa_sites=std_pa_sites,
              pa_ws=pa_ws,    # weights for PA sites
              noise_ws=noise_ws, # weights for noise component
              
              # parameters for fragment length
              mu_frag = mu_frag, # mean of fragment size
              std_frag = std_frag, # standard deviation of fragment size
              dispersion_rate = dispersion_rate, # dispersion_rate x std_frag is actual std
              
              # parameters for read 1 length
              mu_read_utr = mu_read_utr, # mean of r1 length on UTR
              std_read_utr = std_read_utr, # std of r1 length on UTR
              
              # parameters for read 2 length
              mu_read_polya = mu_read_polya, # mean of r2 length on polya
              std_read_polya = std_read_polya, # std of r2 length of polya
              map_rate = map_rate, # probability of r2 that can be mapped to UTR
              
              n_frag=n_frag))
  
}

sim_plot = function(res) {
  ############  plot simulated data ######################
  
  L = res$L
  LA = res$LA
  x = res$x
  l = res$l
  r = res$r
  s = res$polya_part_len_arr
  pa_site_arr = res$theta_arr
  ws = res$ws
  n_frag = res$n_frag
  mu_pa_sites = res$mu_pa_sites
  std_pa_sites = res$std_pa_sites
  
  coverage_cnt = rep(0, L + LA)
  for (i in seq(n_frag)) {
    #cat("i=",i,"\n")
    coverage_cnt[x[i]:(x[i] + l[i] - 1)] = coverage_cnt[x[i]:(x[i] + l[i] -
                                                                1)] + 1
    if (is.na(r[i])) {
      next
    }
    coverage_cnt[(L + s[i] - r[i] + 1):(L + s[i])] = coverage_cnt[(L + s[i] -
                                                                     r[i] + 1):(L + s[i])] + 1
  }
  
  # plot(seq(L+LA),coverage_cnt,type="s",
  #      xlab = paste("UTR | polyA","ws=", paste(sprintf("%.2f",ws), collapse=" ")),
  #      ylab = "coverage count",
  #      main = "red: PA sites, blue: UTR&polyA boarder")
  # abline(v=L, col="blue",lwd=3)
  # for(i in mu_pa_sites){
  #   abline(v=i, col="red",lwd=2,lty=2)
  # }
  
  plot(
    seq(L + LA),
    coverage_cnt,
    type = "s",
    # xlab = paste("UTR | polyA","ws=", paste(sprintf("%.2f",ws), collapse=" ")),
    xlab = paste0("UTR | polyA"),
    ylab = "coverage count",
    main = paste0(paste(
      round(ws, digits = 2), sep = " ", collapse = " "
    )) ,
    xaxt = 'n'
  )
  abline(v = L, col = "#ffed6f", lwd = 3)
  
  ymax = par("usr")[4]
  
  for (i in seq(length(mu_pa_sites))) {
    st = mu_pa_sites[i] - std_pa_sites[i]
    en = mu_pa_sites[i] + std_pa_sites[i]
    x = c(st, st, en, en)
    y = c(0, ymax, ymax, 0)
    lines(
      x,
      y,
      type = 'S',
      col = 'blue',
      lwd = 2,
      lty = 2
    )
    abline(
      v = mu_pa_sites[i],
      col = "green",
      lwd = 2,
      lty = 1
    )
  }
  
  axis(side = 1,
       at = c(seq(0, L + LA, 100)),
       las = 2)
  
  jrt = as.data.frame(table(pa_site_arr[!is.na(pa_site_arr)]))
  loc = as.numeric(levels(jrt[, 1]))
  pmf = jrt[, 2] / sum(jrt[, 2])
  lines(loc,
        ymax * pmf,
        type = 'h',
        col = 'red',
        lwd = 1)
  # legend(
  #   "topleft",
  #   legend = c("pa site (mu)", "pa site (sd)", "pa support", "UTR border"),
  #   col = c("green", "blue", "red", "#ffed6f"),
  #   cex = 0.5,
  #   lty = 1,
  #   box.lty = 0,
  # )
}

# main code
# 
# res = apa_sim(mu_pa_sites = 500,
#               std_pa_sites = 10,
#               pa_ws = 1,
#               map_rate = 0)


import sys
from rpy2.robjects.packages import STAP
import faidx
import gzip
import random
import numpy as np

def main(fasta, bed, outfq, outgt):
  simu = '''run = function(K, alpha_arr){
          source('/home/zhouran/data/proj/2021-0918-apa_evaluation/bin/apamix_simulation.R')
          
          while (TRUE) {
            alpha_list = list(c(9, 1), c(3, 6, 1), c(2, 4, 3, 1))

            ws = MCMCpack::rdirichlet(1, alpha_list[[K]])[1,]
            coverage = round(runif(K) * 5000)
            beta_arr = sample(seq(5, 30, 5), K)
            L = alpha_arr
            # decrease from 100 to 50
            frag_len_mu = 250 + round(50*runif(1))

            frag_len_sd = 20 + round(20*runif(1))

            test <- try({
              tmp_val = apa_sim(
                utr_len = L,
                mu_pa_sites = alpha_arr,
                n_frag = coverage,
                pa_ws = ws[1:K],
                noise_ws = ws[K + 1],
                std_pa_sites = beta_arr,
                mu_frag = frag_len_mu,
                mu_read_utr = 150,
                std_read_utr=0,
                std_frag = frag_len_sd
              )
              tmp_val
            }, silent = TRUE)
            if ((class(test) != "try-error") &&
                (length(test$l) == length(test$x))) {
              break
            }
          }
          return(test)
      }
      '''

  simu = STAP(simu, 'run')
  fa = faidx.Faidx(fasta)


  with open(bed) as fh, \
      gzip.open(outfq,'wb') as r2_out, \
      open(outgt, 'w') as tfh:

      for line in fh:
          line = line.strip().split('\t')
          chrom, st, en, meta, strand = line[:4] + [line[5]]
          simu_dict = simu.run(1, 1000)
          d = { key : list(simu_dict.rx2(key)) for key in simu_dict.names }
          x = d['x']
          l = d['l']
          s = d['s']
          r = d['r']
          beta = d['std_read_utr'][0]
          pa_ws = d['pa_ws'][0]
          noise_ws = d['noise_ws'][0]
          tfh.write(f'{meta}\t{len(x)}\t{pa_ws}\t{noise_ws}\n')

          ig_set = {'nan', 'NA'}
          for inds, s_start in enumerate(x):
              new_header = f'{meta};{x[inds]};{l[inds]};{s[inds]};{r[inds]};{beta}'
              polya = ""
              off_set=0
              if strand == '+':
                  fq_st = int(int(en) - (1000 - s_start + 1))
                  fq_en = int(fq_st + l[inds])

                  if l[inds] + x[inds] > 1000:
                    off_set = abs(int(1000 - (l[inds] + x[inds])))
                    fq_seq = fa.fetch(str(chrom), fq_st + 1, fq_en-off_set, strand)
                    polya = 'A' * off_set
                  else:
                    fq_seq = fa.fetch(str(chrom), fq_st + 1, fq_en, strand)
                  # if not str(s[inds]) in ig_set:
                  #   r1_seq = fa.fetch(str(chrom), int(int(en) - s[inds] + 1), int(en), strand).reverse.complement.seq
                  # else:
                  #   r1_seq = 'N' * 100
              else:
                  fq_st = int(int(st) + (1000 - s_start + 1 ))
                  fq_en = int(fq_st - l[inds])
                  if l[inds] + x[inds] > 1000:
                    off_set = abs(int(1000 - (l[inds] + x[inds])))
                    fq_seq = fa.fetch(str(chrom), fq_en + 1 + off_set, fq_st, strand).reverse.complement
                    polya = 'A' * off_set
                  else:
                    fq_seq = fa.fetch(str(chrom), fq_en + 1, fq_st, strand).reverse.complement

                  # if not str(s[inds]) in ig_set:
                  #   r1_seq = fa.fetch(str(chrom), int(st), int(int(st) + s[inds]), strand).seq
                  # else:
                  #   r1_seq = 'N' * 100

              r2_out.write('''@{}\n{}\n+\n{}\n'''.format(new_header,fq_seq.seq + polya,'I'*(len(fq_seq.seq) + off_set)).encode())
if __name__ == '__main__':
  import fire
  fire.Fire(main)

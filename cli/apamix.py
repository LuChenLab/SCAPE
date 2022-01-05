import os
import sys
import click
import signal
import tqdm
import time
import pandas as pd
import scipy.io, scipy.sparse
from numpy.lib.arraysetops import isin
from pathlib import Path


from loguru import logger
from multiprocessing import Pool
from multiprocessing.pool import MaybeEncodingError
from multiprocessing import set_start_method

from utils.utils import dict_list
from apamix.inference import wraper_process


logger.add('apamix.log',
            rotation='10 MB',
            colorize=True,
            level="DEBUG")


@click.command()
@click.option(
    '--bed',
    type=str,
    help='The target regions (bed format) used for mix inference',
    required=True
    )
@click.option(
    '--bam',
    type=str,
    help='The bam file (sorted and indexed)',
    required=True
    )
@click.option(
    '--out',
    '-o',
    type=str,
    help='The output path',
    required=True
    )
@click.option(
    '--cores',
    type=int,
    help='Num (cores) of region are infering at once',
    default=1
    )
@click.option(
    '--cb',
    type=str,
    help='The cell barcode file, one cb for one line.',
    required=True
    )
@click.option(
    '--tag',
    type=str,
    default='CB,UB',
    help='The cell barcode and UMI tag, for 10X: CB,UB.',
    required=True
    )

@click.option(
    '--n_max_apa',
    type=int,
    default='5',
    help='The maximum number of pA sites. Default value is 5.',
    required=True
    )

@click.option(
    '--n_min_apa',
    type=int,
    default='1',
    help='The minimum number of pA sites. Default value is 1.',
    required=True
    )

@click.option(
    '--max_utr_len',
    type=int,
    default=6000,
    help='The maximum length of UTR. Default value is 6000.'
    )

@click.option(
    '--la_dis_arr',
    type=str,
    default=None,
    help='The distinct lengths of polyA lengths in the dataset. Default: np.arange(self.min_LA, self.max_LA, 10). User counld pass `[10, 30, 50, 70, 90, 110, 130]`'
    )

@click.option(
    '--pmf_la_dis_arr',
    type=str,
    default=None,
    help='The the number of reads for each distinct polyA length. .Default: Unif(min_LA, max_LA). [309912, 4107929, 802856, 518229, 188316, 263208, 101]'
    )


@click.option(
    '--verbose',
    '-v',
    is_flag=True,
    help='Verbose mode'
)
def apamix(
    bed,
    bam,
    out,
    cb,
    tag,
    cores,
    n_max_apa,
    n_min_apa,
    max_utr_len,
    la_dis_arr,
    pmf_la_dis_arr,
    verbose
    ):
    if not all([bed, bam, out, cb]):
        cli(['apamix', '--help'])
        sys.exit(1)

    if not os.path.exists(os.path.join(out, 'tmp')):
        os.makedirs(os.path.join(out, 'tmp'))
        os.makedirs(os.path.join(out, 'huge'))
        os.makedirs(os.path.join(out, 'TooLongRegion'))
        os.makedirs(os.path.join(out, 'TimeConsulting'))

    target_region = open(bed, 'r')
    res_lst = []
    cb_df = pd.read_csv(cb, names=['cb'])
    # cb_df.cb = list(map(lambda x : x.split('-')[0], cb_df.cb.values))

    if la_dis_arr and pmf_la_dis_arr:
        la_dis_arr = eval(la_dis_arr)
        pmf_la_dis_arr = eval(pmf_la_dis_arr)
        if not (isinstance(la_dis_arr, list) and isinstance(pmf_la_dis_arr, list)):
            sys.exit('la_dis_arr and pmf_la_dis_arr were not list format in python.')
    else:
        la_dis_arr = None
        pmf_la_dis_arr = None

    peak_lst = []

    for i in target_region:
        if not i:
            continue
        chrom, st, en, strand = i.strip().split('\t')
        if int(en) - int(st) + 1 > max_utr_len:
            logger.info(f'Skip more than {max_utr_len} UTR, {chrom}_{st}_{en}_{strand}')
            continue

        peak_lst.append(i.strip())

    target_region.close()

    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)

    pool = Pool(processes=cores)
    signal.signal(signal.SIGINT, original_sigint_handler)

    res_lst = []
    try:
        for x in range(len(peak_lst)):
            arg = [peak_lst[x], bam, cb_df, out, tag, verbose, n_max_apa, n_min_apa, la_dis_arr, pmf_la_dis_arr]
            res_lst.append(pool.apply_async(wraper_process, (arg,)))

    except KeyboardInterrupt:
        logger.info('Caught KeyboardInterrupt, terminating workers')
        pool.terminate()
        sys.exit('KeyboardInterrupt')

    pool.close()
    pool.join()

    logger.info('Concating your final sheet')

    md, hd='w', True
    for df in res_lst:
        try:
            df = df.get()
        except MaybeEncodingError:
            continue

        if not isinstance(df, pd.DataFrame):
            continue

        df = df.transpose()
        if hd:
            df.columns = cb_df.cb.values.tolist()
        df.to_csv(f'{out}/pasite.csv.gz', mode=md, header=hd, compression='gzip')
        md, hd='a', False
    logger.info('All done')

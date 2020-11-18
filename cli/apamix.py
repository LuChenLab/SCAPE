import os
import sys
import click
import signal
from numpy.lib.arraysetops import isin
import tqdm
import time
from pathlib import Path
import pandas as pd

from loguru import logger
import scipy.io, scipy.sparse
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
    cores,
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

    peak_lst = []

    for i in target_region:
        if not i:
            continue
        chrom, st, en, strand = i.strip().split('\t')
        if int(en) - int(st) + 1 > 6000:
            too_long_file = f'{out}/TooLongRegion/{chrom}_{st}_{en}_{strand}'
            Path(too_long_file).touch()
            continue

        peak_lst.append(i.strip())

    target_region.close()

    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)

    pool = Pool(processes=cores)
    signal.signal(signal.SIGINT, original_sigint_handler)

    res_lst = []
    try:
        for x in range(len(peak_lst)):
            arg = [peak_lst[x], bam, cb_df, out, verbose]
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

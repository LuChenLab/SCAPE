import sys
import pysam
import fire
from collections import defaultdict
import string, random



def main(bamfile):
    bam_fh = pysam.AlignmentFile(bamfile, 'rb', check_sq=False)
    chars=string.digits + string.ascii_letters + string.punctuation

    outbamfile = pysam.AlignmentFile(bamfile.replace('.bam','.addtag.bam'), 'wb', template = bam_fh)
    add_tag = 'A' * 12
    n = 0
    umi_tag = set()

    for read in bam_fh.fetch(until_eof=True):
        n += 1
        read.set_tag("CB", add_tag, value_type = "Z", replace=True)
        c_umi = ''.join(random.sample(chars, 16))
        while c_umi in umi_tag:
            c_umi = ''.join(random.sample(chars, 16))
        read.set_tag("UR", c_umi, value_type = "Z", replace=True)

        outbamfile.write(read)
        if n % 10000 == 0:
            print(f"processed {n} lines.")
    outbamfile.close()
    bam_fh.close()

if __name__ == '__main__':
 fire.Fire(main)
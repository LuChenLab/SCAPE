import sys
import fire
import pysam
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter, defaultdict

def main(bamfile, config, outprefix):
	bamfile = pysam.AlignmentFile(bamfile,'rb')
	fh = open(config,'r').read().split('\n')
	oh = open(f"{outprefix}.txt", "w")
	align_gap_list = []
	for line in fh:
		if not line: continue

		line = line.split('\t')
		chrom, left_site, right_site, strand = line[:4]

		left_site = int(left_site)
		right_site = int(right_site)
		res = defaultdict(lambda :defaultdict(list))


		for line in bamfile.fetch(chrom, int(left_site) - 1, int(right_site) + 1):
			cigars = set(map(lambda x: x[0], line.cigar))
			if 3 in cigars:
				continue

			start_site = line.reference_start + 1
			end_site = line.reference_end + 1
			query_name = line.query_name

			if line.is_read1:
				res[query_name]['r1'].append(line)
			else:
				res[query_name]['r2'].append(line)
		
		for k,v in res.items():
			label = '_'.join(sorted(v.keys()))
			if label != "r1_r2":
				continue
			r1, r2 = v['r1'], v['r2']
			if len(r1) != 1 or len(r2) != 1:
					continue
			r1, r2 = r1[0], r2[0]
			if strand == "+":
				if r2.is_reverse:
					continue

				discard_len = r1.query_length - r1.query_alignment_length + r2.query_length - r2.query_alignment_length
				align_gap = r1.reference_end - r2.reference_start + 1 + len(r1.get_tag("DT")) + discard_len
				if align_gap < 0:
					continue
				align_gap += discard_len
				
			else:
				"""
				==r1==>
				   <==r2==
				"""
				if not r2.is_reverse:
					continue

				discard_len = r1.query_length - r1.query_alignment_length + r2.query_length - r2.query_alignment_length
				align_gap = r1.reference_start - r2.reference_end + 1 + len(r1.get_tag("DT"))
				if align_gap < 0:
					continue
				align_gap += discard_len

			if align_gap > 1000:
				continue
			oh.write("{}".format(align_gap) + "\n")
			align_gap_list.append(align_gap)

	oh.close()
	points = sns.kdeplot(np.array(align_gap_list)).get_lines()[0].get_data()
	x, y = points[:2]
	plt.axvline(300, color='r')
	plt.savefig(f"{outprefix}.pdf")

			

if __name__ == '__main__':
	fire.Fire(main)

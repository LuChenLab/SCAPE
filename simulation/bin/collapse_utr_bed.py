import sys
from collections import defaultdict

def main(args):
	bedin, bedout = args
	out_dic = defaultdict(list)
	with open(bedin) as fh, open(bedout, 'w') as fo:
		for line in fh:
			gid = line.strip().split('\t')[-1]
			out_dic[gid].extend([line])

		for gid, detail in out_dic.items():

			detail = list(map(lambda x: x.strip().split('\t'), detail))
			detail_ = list(map(lambda x: ','.join(x), zip(*detail)))
			chrom, st, en = detail_[:3]
			strand = detail_[-2].split(',')[0]
			chrom = chrom.split(',')[0]
			if strand == '+' or strand == '-':
				st = str(min(map(lambda x: int(x), st.split(','))))
				en = str(max(map(lambda x: int(x), en.split(','))))
			else:
				raise(f'Unrecognize strand {strand}')
			for i in detail:
				fo.write('\t'.join([chrom, st, en, strand] + i) + '\n')


if __name__ == '__main__':
	main(sys.argv[1:])

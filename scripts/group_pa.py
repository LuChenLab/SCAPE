import sys
import pandas as pd

def make_df(file, label):
	df = pd.read_table(file, header=None, names = ['pa_site'])
	df['label'] = label
	tmp_df = pd.DataFrame(list(map(lambda x: x.split(':'), df.pa_site.values.tolist())), columns = ['chrom', 'pa','beta','strand'])
	df = pd.concat([df, tmp_df],axis=1)
	return(df)

def main(args):
	files, labels, outfile, group_threshold = args
	group_threshold = int(group_threshold)
	tmp_list = []
	for file, label in zip(files.split(','), labels.split(',')):
		tmp_list.append(make_df(file, label))
	df = pd.concat(tmp_list)
	df.pa = df.pa.astype('int64')
	df = df.sort_values(['chrom','pa']).groupby(['chrom', 'strand'])
	res_lst = []
	ind = 0
	for group_info, group_df in df:
		chrom, strand = group_info
		group_ids = (group_df.pa > group_df.pa.shift() + group_threshold).cumsum()
		group_df = group_df.groupby(group_ids)
		for _, sub_df in group_df:
			ind += 1
			# print(ind)
			if ind % 10000 == 0:
				print(f'prcessing {ind} line')
			pasite = sub_df.pa.values.tolist()
			sub_df['collapse_pa'] = pasite[0]  if strand == '-' else pasite[-1]
			res_lst.append(sub_df)
	df = pd.concat(res_lst)
	df.to_csv(outfile, sep = '\t', index=False)

if __name__ == '__main__':
	main(sys.argv[1:])
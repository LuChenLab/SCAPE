import sys
import click
import pandas as pd
import numpy as np



def make_df(file, label):
	df = pd.read_csv(file)
	df = df.rename(columns = {df.columns[0]:'pa_site'})
	df.loc[:,'pseudo_expr'] = df.drop('pa_site', axis=1).sum(axis=1)
	df = df[['pa_site', 'pseudo_expr']]
	df.loc[:,'label'] = label
	tmp_df = pd.DataFrame(list(map(lambda x: x.split(':'), df.pa_site.values.tolist())), columns = ['chrom', 'pa','beta','strand'])
	df = pd.concat([df, tmp_df],axis=1)
	return(df)


def run(files, labels, outfile, group_threshold):
	group_threshold = int(group_threshold)
	tmp_list = []
	for file, label in zip(files.split(','), labels.split(',')):
		tmp_list.append(make_df(file, label))
	df = pd.concat(tmp_list)
	df.pa = df.pa.astype('int64')
	df = df.sort_values(['chrom', 'pa']).groupby(['chrom', 'strand'])
	res_lst = []
	ind = 0
	for group_info, group_df in df:
		chrom, strand = group_info
		group_ids = (group_df.pa > group_df.pa.shift() + group_threshold).cumsum()
		group_df = group_df.groupby(group_ids)
		for _, sub_df in group_df:
			ind += 1
			if ind % 10000 == 0:
				print(f'Prcessed {ind} line')

			#  skip the SettingWithCopyError information
			sub_df = sub_df.copy()
			pasite = sub_df.pa.values.tolist()
			pseudo_exp = sub_df.pseudo_expr.values
			if sum(pseudo_exp) != 0:
				sub_df.loc[:, 'collapse_pa'] = pasite[np.argmax(sub_df.pseudo_expr.values)]
			else:
				sub_df.loc[:, 'collapse_pa'] = pasite[0] if strand == '-' else pasite[-1]
			res_lst.append(sub_df)

	df = pd.concat(res_lst)
	df.to_csv(outfile, sep = '\t', index=False)

@click.command()
@click.option('--files',
			  help='Input files to group pA, e.g. sample1/pasite.csv.gz,sample2/pasite.csv.gz',
			  type=str
			  )
@click.option('--labels',
			  help='Sample label for each input file, sample1,sample2',
			  type=str
			  )
@click.option('--outfile',
			  help='output file')
@click.option('--group_threshold',
			  help='The minimum distance for grouping pA, defualt: 100',
			  type=int,
			  default=100
			  )
def cli(files, labels, outfile, group_threshold):
	run(files,
		labels,
		outfile,
		group_threshold
		)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		cli(['--help'])
		sys.exit(0)
	cli()

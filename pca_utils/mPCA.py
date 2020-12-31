#!/usr/bin/env python

import sys
import argparse as ap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
from sklearn import decomposition
from sklearn import preprocessing as prep
from sklearn.metrics import pairwise_distances
from sklearn import manifold


def read_dists(fin, perc=98):
	df = pd.read_table(fin, header=0, index_col=0)
	df = df / df.max()
	p = np.percentile(df, perc)
	df[df > p] = p
	df = df / df.max()

	return df


def pca(f):
	pca = decomposition.PCA(n_components=2)
	pca = pca.fit(f)
	ft = pca.transform(f)

	return (ft, pca.explained_variance_ratio_)


def pcoa(f, metric):
	dists = pairwise_distances(f, metric=metric)
	mds = manifold.MDS(n_components=2, max_iter=5000, eps=1e-9, dissimilarity='precomputed')
	pos = mds.fit(dists).embedding_

	pos *= np.sqrt((f**2).sum()) / np.sqrt((pos**2).sum()) #rescale the data

	#rotate the data
	pca = decomposition.PCA(n_components=2)
	pca = pca.fit(pos)
	ft = pca.transform(pos)

	return (ft, pca.explained_variance_ratio_)


def draw(ft, l, m, pca_var, fout):
	l = l.values.flatten().astype('float')
	m = m.values.flatten().astype('str')
	fig, ax = plt.subplots()

	for i in np.unique(m):
		ax.scatter(ft[m==i, 0], ft[m==i, 1], c=l[m==i], marker=i, s=200, cmap=plt.cm.jet, vmin=min(l), vmax=max(l), edgecolors='None', alpha=0.6)

	ax.set_xlabel('PC1 (' + str(np.round(100*pca_var[0], decimals=2)) + '%)')
	ax.set_ylabel('PC2 (' + str(np.round(100*pca_var[1], decimals=2)) + '%)')
	ax.set_title('PCoA')
	ax.set_xlim([min(ft[:, 0]), max(ft[:, 0])])
	ax.set_ylim([min(ft[:, 1]), max(ft[:, 1])])
	ax.legend()
	fig.savefig(fout, bbox_inches='tight')


def read_params(args):
	parser = ap.ArgumentParser(description='Microbiome Principal Component Analysis')
	arg = parser.add_argument
	arg('inp_f', metavar='INPUT_FILE', nargs='?', default=sys.stdin, type=str, help="the input dataset file [stdin if not present]")
	arg('out_fig', metavar='OUT_FIGURE', nargs='?', default=None, type=str, help="the output figure file")
	arg('-z', '--feature_identifier', type=str, default='k__', help="the feature identifier")
	arg('-d', '--define', type=str, help="define discrete labels")
	arg('-d2', '--define2', type=str, help="define continuous labels")
	arg('-dm', '--definem', type=str, help="define discrete markers")
	arg('--precomputed', action='store_true', default=False, help="whether the input file is a distance matrix or not (default false)")
	arg('--pcoa', action='store_true', default=False, help="plot PCoA (default PCA)")
	arg('-m', '--metric', default='euclidean', help="define distance metric for PCoA")

	par = vars(parser.parse_args())

	return par


if __name__ == "__main__":
	par = read_params(sys.argv)

	if par['precomputed']:
		ft = read_dists(par['inp_f'])
	else:
		f = pd.read_csv(par['inp_f'], sep='\t', header=None, index_col=0, dtype=unicode)
		f = f.T

		if par['define']:
			d = pd.DataFrame([s.split(':') for s in par['define'].split(',')])
			l = pd.DataFrame([0]*len(f))
			for i in range(len(d)):
				l[(f[d.iloc[i,1]].isin(d.iloc[i,2:])).tolist()] = d.iloc[i,0]
		elif par['define2']:
			l = pd.DataFrame(f[par['define2']])
			l = l.convert_objects(convert_numeric=True).fillna(value=0.0)
		else:
			le = prep.LabelEncoder()
			le.fit(f.iloc[:,0])
			l = pd.DataFrame(le.transform(f.iloc[:,0]))

		m = pd.DataFrame(['o']*len(f))

		if par['definem']:
			d = pd.DataFrame([s.split(':') for s in par['definem'].split(',')])
			for i in range(len(d)):
				m[(f[d.iloc[i,1]].isin(d.iloc[i,2:])).tolist()] = d.iloc[i,0]

		feat = [s for s in f.columns if sum([s2 in s for s2 in par['feature_identifier'].split(':')]) > 0]
		f = f.loc[:, feat].astype('float')

		if par['pcoa']:
			ft, pca_var = pcoa(f.values, par['metric'])
		else:
			ft, pca_var = pca(f.values)

	assert ft.all()
	draw(ft, l, m, pca_var, par['out_fig'])

#!/usr/bin/env python3
# Author: Bhuvan Molparia

import os
import argparse

import numpy as np

from parse_tax_summary import tax_parser

'''
This script calculates the alpha diversity and the beta diversity for
any given taxonomic distribution at TAX LEVEL 3
Inputs are - Taxonomy file (output of Mothur), sample sheet, outfile path
'''

def parse_sample(filepath):
	fin  = open(filepath,'r')
	data = fin.readlines()
	fin.close()

	snames = []
	sdict  = {}
	for line in data:
		line = line.replace('\n','').split('\t')
		sdict[line[1]] = line[0]
		snames.append(line[1])

	return snames,sdict

class MethodError(Exception):
	pass

def calc_alpha(counts,method='shannon'):

	counts = list(filter(lambda a: a != 0, counts))
	fracs  = np.array(counts,dtype=float)
	fracs  = fracs/sum(fracs)

	if method == 'shannon':
		alpha = -sum(fracs*np.log(fracs))
	elif method == 'simpson':
		alpha = sum(fracs*fracs)
	else:
		raise MethodError('Unknown method: '+method)

	return alpha

def calc_beta(counts1,counts2,method='fractions'):

	counts1 = np.array(counts1,dtype=float)
	counts2 = np.array(counts2,dtype=float)

	if method == 'fractions':

		fracs1 = counts1/sum(counts1)
		fracs2 = counts2/sum(counts2)
		BCdiss = sum(np.absolute(fracs1-fracs2))/2

	elif method == 'counts':
		BCdiss = sum(np.absolute(counts1-counts2))/(sum(counts1)+sum(counts2))
	else:
		raise MethodError('Unknown method: '+method)

	return BCdiss

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='This script calculates a '
	'variety of species richness metrics for a given taxonomy summary file '
	'produced by mothur')
	parser.add_argument('-i', dest='inp', help='path to the input data file')
	parser.add_argument('-s', dest='sam', help='path to optional tsv file '
	'conataining mappings of column names in input data to sample identifiers',
	default = None)
	parser.add_argument('-o', dest='out', help='prefix to files for output '
	' files')

	args = parser.parse_args()

	data  = tax_parser(args.inp)

	if args.sam == None:
		sample_names = sorted(list(data))
		sample_dict = {key:key for key in sample_names}
	else:
		sample_names, sample_dict = parse_sample(args.sam)

	## ------------- Alpha Index ------------------------
	outf1_path = '{}.alpha_diversity.txt'.format(args.out)
	with open(outf1_path,'w') as outf1:
		outf1.write('Sample\tShannonInd\tSimpsonInd\n')
		for i in sample_names:
			samI = sample_dict[i]
			counts  = data[samI][3]['counts']
			shannon = str(calc_alpha(counts))
			simpson = str(calc_alpha(counts,method='simpson'))
			outf1.write(i+'\t'+shannon+'\t'+simpson+'\n')

	## ------------- Beta Index -------------------------
	outf2_path = '{}.beta_diversity.txt'.format(args.out)
	with open(outf2_path,'w') as outf2:

		outf2.write('Bray-Curtis Dissimilarity - Fractions\n\n')
		outf2.write('\t'+('\t').join(sample_names)+'\n')
		for i in sample_names:
			samI = sample_dict[i]
			outf2.write(i+'\t')
			counts1 = data[samI][3]['counts']

			for j in sample_names:
				samJ = sample_dict[j]
				counts2 = data[samJ][3]['counts']
				BCDiss  = str(calc_beta(counts1,counts2))
				outf2.write(BCDiss+'\t')

			outf2.write('\n')

		outf2.write('\n\n Bray-Curtis Dissimilarity - Counts\n\n')
		outf2.write('\t'+('\t').join(sample_names)+'\n')
		for i in sample_names:
			samI = sample_dict[i]
			outf2.write(i+'\t')
			counts1 = data[samI][3]['counts']

			for j in sample_names:
				samJ = sample_dict[j]
				counts2 = data[samJ][3]['counts']
				BCDiss  = str(calc_beta(counts1,counts2,method='counts'))
				outf2.write(BCDiss+'\t')

			outf2.write('\n')

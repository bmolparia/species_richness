#!/usr/bin/env python3
# Author: Bhuvan Molparia

def tax_parser(fpath):
    ''' Main function to be called to run the parser '''

    fhand = open(fpath,'r')

    header = fhand.readline()
    rank_ind, taxon_ind, total_ind, samples = parse_header(header)

    Data = {}
    for i in samples:
        Data[i] = {}

    dataline = fhand.readline()

    while dataline:
        rank_id, taxon, total, taxlvl, sample_values = parse_line(dataline,
        rank_ind, taxon_ind, total_ind, samples)

        for samp in Data:

            if taxlvl not in Data[samp]:
                Data[samp][taxlvl] = {}
                Data[samp][taxlvl]['taxon']  = []
                Data[samp][taxlvl]['counts'] = []
                Data[samp][taxlvl]['rank_id'] = []

            Data[samp][taxlvl]['taxon'].append(taxon)
            Data[samp][taxlvl]['counts'].append(sample_values[samp])
            Data[samp][taxlvl]['rank_id'].append(rank_id)

        dataline = fhand.readline()

    fhand.close()
    return Data

def hierarchy(rank):
    ''' Function to get the heirarchy in the taxonomy classification'''
    rank = rank.split('.')
    i = 0
    while i < len(rank):
        yield ('.').join(rank[0:i+1])
        i+= 1

def parse_line(line, rank_ind, taxon_ind, total_ind, samples):
    ''' Function to parse indivdual rows and get extract data'''
    line = line.replace('\n', '').split('\t')

    rank_id = line[rank_ind]
    taxon  = line[taxon_ind]
    total  = int(line[total_ind])
    taxlvl = rank_id.count('.')

    sample_values = {}
    for i in samples:
        sample_values[i] = int(line[samples[i]])

    return rank_id, taxon, total, taxlvl, sample_values

def parse_header(line):
    ''' Function to parse the header of a mothur taxonomy summary output file
    Sample header:
    taxlevel   rankID  taxon   daughterlevels  total  v34f0r0 v34f10r10 ...
    '''
    line = line.replace(' ', '').replace('\n', '').split('\t')

    rank_ind  = line.index('rankID')
    taxon_ind = line.index('taxon')
    total_ind = line.index('total')

    samples = {}
    for i in line[total_ind+1:]:
        if i != '':
            samples[i]  = line.index(i)

    return rank_ind, taxon_ind, total_ind, samples


if __name__ == "__main__":

    import sys
    fpath = sys.argv[1]

    d = tax_parser(fpath)
    for i in d:
        print(i, '\t', d[i])

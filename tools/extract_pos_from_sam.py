#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 17:48:12 2016

@author: antlaplante
"""

import argparse
import pandas as pd
from Bio import SeqIO


def get_aln(header):

    (name, _, chrom, pos, qual) = header

    return {'name' : name, \
            'chrm' : chrom, \
            'pos' : int(pos), \
            'qual' : int(qual)}

def read_id2idx(reads_fn):

    if reads_fn[-1] == 'q':
        fmt = 'fastq'
    else:
        fmt = 'fasta'

    reads_idx={}
    reads_fh = open(reads_fn, "rU")
    num_read=0
    for record in SeqIO.parse(reads_fh, fmt):
        reads_idx[record.id] = num_read
        num_read += 1
    reads_fh.close()

    return reads_idx

def get_headers(sam_fn):

    # Load headers from sam file
    headers = []
    fh = open(sam_fn, 'rb')
    for line in fh:
        header=line.split('\t')[:5]
        headers.append(header)
    fh.close()

    # Check if eukaryotic or prokaryotic and strip headers
    chr_names = []
    skipline = 0
    for header in headers:
        if header[0] == '@SQ':
            chr_names.append(header[1][3:])
            skipline += 1
        elif header[0][0] == '@':
            skipline += 1
            continue
        else: break

    headers = headers[skipline:]

    return (headers, chr_names)

def algn_dic(headers, reads_id_dic):

    algns = {}

    for header in headers:
        aln = get_aln(header)
        read_idx = reads_id_dic[aln['name']]

        if algns.has_key(read_idx):
            if algns[read_idx]['qual'] > aln['qual']:
                continue

        algns[read_idx] = aln

    return algns

def int_to_roman(input):
    """
    Convert an integer to Roman numerals.
    """

    if type(input) != type(1):
        raise TypeError, "expected integer, got %s" % type(input)
    if not 0 < input < 4000:
        raise ValueError, "Argument must be between 1 and 3999"
    ints = (1000, 900,  500, 400, 100,  90, 50,  40, 10,  9,   5,  4,   1)
    nums = ('M',  'CM', 'D', 'CD','C', 'XC','L','XL','X','IX','V','IV','I')
    result = ""
    for i in range(len(ints)):
        count = int(input / ints[i])
        result += nums[i] * count
        input -= ints[i] * count

    return result


# Inputs
parser = argparse.ArgumentParser()
parser.add_argument("alignments", help="path .sam file from BWA-bwasw or GraphMap")
parser.add_argument("reads", help="path to reads")
parser.add_argument('-o', '--out', default='reads_position.csv', \
                    help='output file with position of reads and chromosome' \
                    'to which they belong')

args = parser.parse_args()

(headers, chr_names) = get_headers(args.alignments)

names_chr = {name : int_to_roman(idx+1) for (idx, name) in enumerate(chr_names)}

reads_dic = read_id2idx(args.reads)
algnmts = algn_dic(headers, reads_dic)

for (idx, read_aln) in algnmts.items():
    read_aln['number'] = idx
    if names_chr.has_key(read_aln['chrm']):
        read_aln['chr_nb'] = names_chr[read_aln['chrm']]
    else:
        read_aln['chr_nb'] = '*'

alns_df = pd.DataFrame(algnmts.values()).sort_values('number')
alns_df.to_csv(args.out, sep='\t')

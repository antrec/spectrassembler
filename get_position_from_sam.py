#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 17:48:12 2016

@author: antlaplante
"""


import sys
import argparse
import csv
from Bio import SeqIO
import numpy as np

# define alignment information class
class Aln:
    """ Alignment class
    input : BWA bwasw alignment header """
    def __init__(self, header):
        read_name, flag, chromosome, position, quality = header
        self.name = read_name
        self.chrm = chromosome
        self.pos = int(position)
        self.qual = int(quality)
# Inputs
parser = argparse.ArgumentParser()
parser.add_argument("samfile", help="path .sam file from BWA-bwasw")
parser.add_argument("fastafile", help="path to reads in fasta format")
parser.add_argument("-r", "--root", type=str, default="./",
                    help="directory where to store position file.")


args = parser.parse_args()

doChr = False
# sam_fn=sys.argv[1]
# reads_fn=sys.argv[2]

# load sam file headers
headers=[]
fh = open(args.samfile,'r')
for line in fh:
    header=line.split('\t')[:5]
    headers.append(header)
fh.close()

num_chr=0
chr_names=[]
for header in headers:
    if header[0]=='@SQ':
        num_chr += 1
        chr_names.append(header[1][3:])
    else: break

if num_chr > 1:
    doChr = True
headers=headers[num_chr:]

# make dictionary read name --> read number
reads_idx={}
ref=open(args.fastafile, "rU")
num_read=0
for record in SeqIO.parse(ref,"fasta"):
    reads_idx[record.id] = num_read
    num_read += 1
ref.close()


# keep best alignment for each read aligned
aln_prec = Aln(headers[0])
best_alns={}
for header in headers:
    aln = Aln(header)
    if not(aln.name == aln_prec.name and aln.qual < aln_prec.qual):
        best_alns[aln.name] = aln
        aln_prec = aln

# write to file(s)
if doChr:
    chr_fn ='chr_clusters.csv'
    chr_fh = open(chr_fn,'wb')
    chr_writer = csv.writer(chr_fh)
    chr_num = 0
for chr_name in chr_names:
    if doChr:
        chr_num += 1
        reads_in_chr = [read_name for read_name in best_alns.keys() if best_alns[read_name].chrm == chr_name]
        reads_in_chr_idx = np.array([reads_idx[read_name] for read_name in reads_in_chr])
        positions = np.array([best_alns[read_name].pos for read_name in reads_in_chr])
        order = np.argsort(reads_in_chr_idx)
        chr_writer.writerow(reads_in_chr_idx[order]+1)
        pos_fn = args.root + '/euk_chr_%d.csv' % (chr_num)
    else:
        reads_in_chr = best_alns.keys()
        reads_in_chr_idx = np.array([reads_idx[read_name] for read_name in reads_in_chr])
        positions = np.array([best_alns[read_name].pos for read_name in reads_in_chr])
        order = np.argsort(reads_in_chr_idx)
        pos_fn = args.root + '/reads_position.csv'

    pos_fh = open(pos_fn,'wb')
    pos_writer = csv.writer(pos_fh)
    pos_writer.writerow(positions[order])
    pos_fh.close()

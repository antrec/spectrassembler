#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 11:58:50 2015

@author: Antoine Recanati at INRIA.
Adapted from nanocorrect [https://github.com/jts/nanocorrect/] (Loman, Quick, Simpson), (c) 2015 Ontario Institute for Cancer Research.

INPUT: consensus in sliding windows from POA
OUTPUT : whole consensus in the contig

"""
import argparse
import subprocess
import os
from Bio import AlignIO

def clustal2consensus(fn):

    alignment = AlignIO.read(fn, "clustal")
    consensus_row = -1

    for (i, record) in enumerate(alignment):
        if record.id == 'CONSENS0':
            consensus_row = i

    if consensus_row == -1:
        return ""

    # Extract the consensus sequence
    consensus = str(alignment[consensus_row].seq)
    consensus = consensus.replace('-', '')

    return consensus

def run_poa_and_consensus(in_fn, out_fn, poa_path, poa_mat_path):

    if poa_path == None:
        poa_path = 'poa'
    if poa_mat_path == None:
        poa_mat_path = 'poa-score.mat'
    elif poa_mat_path[-13:] != 'poa-score.mat':
        poa_mat_path += '/poa-score.mat'

    cmd = "%s -read_fasta %s -clustal %s -hb %s" % (poa_path, in_fn, out_fn, poa_mat_path)
    p = subprocess.Popen(cmd, shell=True)
    p.wait()
    consensus =  clustal2consensus(out_fn)
    os.remove(in_fn)
    os.remove(out_fn)

    return consensus

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--root", type=str, default="./",
    help="directory where to look for files \
(the one given to spectral_layout_from_minimap.py).")
parser.add_argument("--poa_path",
    help="path/to/poa/poa (path to the executable)")
parser.add_argument("--poa_mat_path", help="path/to/poa-score.mat")
parser.add_argument("-cc", "--contig", type=str, required=True,
    help="index of contig you wish to compute consensus for")
parser.add_argument("--w_len", type=int, default=2500,
    help="length of consensus windows for POA.")
parser.add_argument("--w_ovl_len", type=int, default=1250,
    help="overlap length between two successive \
consensus windows.")
parser.add_argument("-v", "--verbosity", action="count", default=0,
    help="verbosity level (nothing : 0, -v : 1; -vv : 2)")

# Assign variables and constants
args = parser.parse_args()
root = args.root
cc_idx = args.contig
w_len = args.w_len
w_ovl_len = args.w_ovl_len
poa_path = args.poa_path
poa_mat_path = args.poa_mat_path
verb = args.verbosity

# Count number of windows
try:
    cmd = "ls %s/cc_%d/poa_in_cc_%d_win_*.fasta.clustal.out" % (root, cc_idx, cc_idx)
    n_win = len(subprocess.check_output(cmd, shell=True))
except:
    n_win = 10000 # quick fix in case of problem with output of subprocess

# Write the consensus sequences from all windows to a file
poa_in_fn = "%s/poa_in_cons_cc_%d.fasta" % (root, cc_idx)

# Initialize
fn = "%s/cc_%d/poa_in_cc_%d_win_%d.fasta.clustal.out" % (root, cc_idx, cc_idx, 0)
win_cons_seq = clustal2consensus(fn)
whole_cons = win_cons_seq

# Incrementally add consensus between window k and window k+1
for w_idx in xrange(1, n_win):
    try:
        fn = "%s/cc_%d/poa_in_cc_%d_win_%d.fasta.clustal.out" % (root, cc_idx, cc_idx, w_idx)
        next_win_seq = clustal2consensus(fn)
        next_win_len = len(next_win_seq)
        whole_cons_len = len(whole_cons)
        kept_len = max(0, whole_cons_len - next_win_len - margin)
        cons0 = whole_cons[:kept_len]
        cons1 = whole_cons[kept_len:]

        # Write end of current consensus long sequence and next consensuns \
#window sequence in poa_in file
        poa_in_fh = open(poa_in_fn, "wb")
        poa_in_fh.write(">end_of_current_cons\n%s\n" % (cons1))
        poa_in_fh.write(">cons_in_window_%d\n%s\n" % (w_idx, next_win_seq))
        poa_in_fh.close()

        #Run poa to include next
        out_fn = "%s/cc_%d/poa_out_cons_cc%d_win_%d" % (root, cc_idx, cc_idx, w_idx)
        cons1b = run_poa_and_consensus(poa_in_fn, out_fn)
        whole_cons = cons0 + cons1b

    except:
        if verb >= 2:
            print("Error for file %s" % (fn))

if verb >= 1:
    print("extracted and merged sequences in windows for contig %d. Consensus length %dbp" % \
     (cc_idx, len(whole_cons)))

consensus_fn = "consensus_cc_%d.fasta" % (cc_idx)
consensus_fh = open(consensus_fn, "wb")
consensus_fh.write(">consensus_from_windows_contig_%d\n%s\n" % (cc_idx, whole_cons))
consensus_fh.close()

# Clean up
os.remove(poa_in_fn)

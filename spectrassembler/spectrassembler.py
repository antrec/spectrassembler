#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""
Spectrassembler main program

@author: Antoine Recanati
"""
from __future__ import print_function
from time import time
import sys
import argparse
from functools import partial
from multiprocessing import Pool
import numpy as np
from Bio import SeqIO
from scipy.sparse import coo_matrix
from scipy.stats.mstats import mquantiles

from overlaps import compute_positions, compute_overlaps
from spectral import sym_max, remove_bridge_reads, reorder_mat_par, reorder_mat
from consensus import run_spoa_in_cc, merge_windows_in_cc
from ioandplots import fill_args_opts, make_dir, oprint, write_layout_to_file, plot_cc_pos_v_ref

# Parse arguments and define global variables
t0 = time()

parser = argparse.ArgumentParser(description="De novo experimental assembler"
                                             "based on a spectral algorithm to reorder the reads")
parser.add_argument("-r", "--root", default="./",
                    help="directory where to store layout and consensus files.")
parser.add_argument("-f", "--READS_FN", required=True,
                    help="path to reads file (fasta or fastq)")
parser.add_argument("-m", "--minimapfn", required=True,
                    help="overlap file path (from minimap in PAF format).")
parser.add_argument("--min_cc_len", type=int, default=10,
                    help="minimum number of reads for a contig to be considered")
parser.add_argument("--w_len", type=int, default=3000,
                    help="length of consensus windows for POA.")
parser.add_argument("--w_ovl_len", type=int, default=2000,
                    help="overlap length between two successive consensus windows.")
parser.add_argument("--len_thr", type=int, default=3500,
                    help="threshold on length of overlaps (similarity matrix preprocessing).")
parser.add_argument("--sim_qtile", type=float, default=0.4,
                    help="quantile threshold on overlap score (similarity matrix preprocessing.)" \
                         "0.5 means you keep only overlaps with num_match > quantile(num_matches, 0.5)")
parser.add_argument("-v", "--verbosity", action="count", default=1,
                    help="verbosity level (-v, -vv or none)")
parser.add_argument("--ref_pos_csvf",
                    help="csv file with position of reads (in same order as in READS_FN)" \
                         "obtained from BWA, in order to plot reads position found vs reference.")
parser.add_argument("--spoapath", default="tools/spoa/spoa",
help="path to spoa executable")
parser.add_argument("--nproc", help="number of parallel processes", type=int,
                    default=1)
parser.add_argument("--margin", type=int, default=1250,
                    help="number of bases to add to current consensus to make sure it overlaps next window")
parser.add_argument("--trim_margin", type=int, default=200,
                    help="length to cut in beginning and end of consensus sequences from spoa (where the consensus is" \
                         "less good)")
parser.add_argument("--julia", default=None,
help="path to Julia (optional,"\
"though eigenvector computations are clearly faster in Julia than in Python)")


args = parser.parse_args()
opts = fill_args_opts(args)
ROOT_DIR = opts['ROOT_DIR']
VERB = opts['VERB']

# Load reads
reads_fh = open(args.READS_FN, "rU")
record_list = list(SeqIO.parse(reads_fh, opts['READS_FMT']))
reads_fh.close()
oprint("Reads loaded. Compute overlaps from files...", dt=(time() - t0), cond=(VERB >= 2))

# Compute overlaps from the files
(read_nb2id, ovl_list, I, J, K, num_match, ovl_len, n_reads) = compute_overlaps(args.minimapfn, record_list)

# Threshold based on overlaps value (number of matches) and length
THR = mquantiles(num_match, args.sim_qtile)
oprint("THR = %1.1f " % THR)
cond1 = (num_match > THR)
cond2 = (ovl_len > opts['LEN_THR'])
idxok = np.argwhere(cond1 * cond2)[:, 0]
num_match_l = num_match
I = I[idxok]
J = J[idxok]
num_match = num_match[idxok]
# ovl_len = ovl_len[idxok]
K = K[idxok]

# Construct similarity matrix
oprint("Construct thresholded similarity matrix...", dt=(time() - t0), cond=(VERB >= 2))
sim_mat = coo_matrix((num_match, (I, J)), shape=(n_reads, n_reads), dtype=int).tocsr()
oprint("Pre-process similarity matrix...", dt=(time() - t0), cond=(VERB >= 2))

# Overlap index array : overlap(i,j) = ovl_list[k], with k = ovl_idx_arr[i,j]
ovl_idx_arr = coo_matrix((K, (I, J)), shape=(n_reads, n_reads), dtype=int).tocsr()
ovl_idx_arr = sym_max(ovl_idx_arr)

# Symmetrize the matrix when it is not already symmetric
sim_mat = sym_max(sim_mat)
# sim_mat = (sim_mat + sim_mat.T)

# Remove "connecting reads"
sim_mat = remove_bridge_reads(sim_mat)
del I, J, K, ovl_len, num_match
oprint("Similarity matrix built and preprocessed. Reorder it with spectral ordering...", dt=(time() - t0),
       cond=(VERB >= 1))

# Reorder connected components with spectral ordering
ccs_list = []
cc = range(sim_mat.shape[0])
qtile = args.sim_qtile
t_start_layout = time()
# reorder_submat(sim_mat, cc, num_match_l, qtile, ccs_list, opts)
thr_list = []
new_qtile = qtile
for k in range(40):
    thr_sub = float(mquantiles(num_match_l, new_qtile))
    thr_list.append(thr_sub)
    new_qtile += min(0.1, 0.5 * (1. - new_qtile))
del num_match_l
if opts['N_PROC'] > 1:
    ccs_list = reorder_mat_par(sim_mat, thr_list, opts)
else:
    ccs_list = reorder_mat(sim_mat, thr_list, opts['MIN_CC_LEN'], opts['VERB'])

t_rough_layout = time() - t_start_layout
oprint("Rough layout computed in %3.3f." % (t_rough_layout),
        dt=(time() - t0), cond=(VERB >= 1))

# Sort by length of connected component
ccs_list.sort(key=len, reverse=True)
oprint("Compute fine grained layout and run spoa in connected components...", dt=(time() - t0), cond=(VERB >= 1))

# If root_dir does not exist, create it
make_dir(ROOT_DIR)

t_total_finegrained = 0

# Get fine-grained layout with dictionary of overlaps in each connected component
for (cc_idx, cc) in enumerate(ccs_list):

    # Restrict overlap index array to reads in the connected component (contig)
    # ovl_idx_cc = ovl_idx_arr.copy().tocsc()[:, cc]
    # ovl_idx_cc = ovl_idx_cc.tocsr()[cc, :]
    ovl_idx_cc = ovl_idx_arr[cc,:][:,cc]
    # symmetrize if the overlapper does not count overlap twice for (i,j) and (j,i)
    # ovl_idx_cc = sym_max(ovl_idx_cc)

    # Compute fine-grained position and strand of each read in connected component
    t_start_fg_layout = time()
    (strand_list, bpos_list, epos_list) = compute_positions(cc, read_nb2id, ovl_list, ovl_idx_cc)
    t_finegrained = time() - t_start_fg_layout
    t_total_finegrained += t_finegrained
    msg = "Positions computed in connected component"\
    "%d/%d in %3.3f.\n Now run spoa if provided." % (cc_idx,
    len(ccs_list) - 1, t_finegrained)
    oprint(msg, dt=(time() - t0), cond=(VERB >= 2))

    # Write file with layout
    layout_fn = "%s/cc%d.layout" % (ROOT_DIR, cc_idx)
    write_layout_to_file(layout_fn, strand_list, bpos_list, epos_list, cc, read_nb2id)
    msg = "layout written to file %s" % (layout_fn)
    oprint(msg, dt=(time() - t0), cond=(VERB >= 2))

    if opts['DO_PLOT_POS_V_REF']:
        msg = "Edit graphic : position of reads found by algorithm vs reference"
        oprint(msg, dt=(time() - t0), cond=(VERB >= 2))
        figpath = ROOT_DIR + "/pos_found_vs_ref_cc%d.eps" % (cc_idx)
        plot_cc_pos_v_ref(opts['REF_POS_CSVF'], cc, bpos_list, figpath)

    # Generate contigs through multiple sequence alignment
    if opts['DO_SPOA']:
    # if False:
        # Compute consensus in windows
        run_spoa_in_cc(record_list, cc_idx, cc, strand_list, bpos_list,
        epos_list, opts)

        if opts['N_PROC'] == 1:
        # Merge windows to get consensus
            cons_in_cc = merge_windows_in_cc(cc_idx, opts)
            print(">contig_%d\n%s" % (cc_idx, cons_in_cc), file=sys.stdout)

            msg = "Consensus computed in connected component %d/%d. " % (cc_idx, len(ccs_list) - 1)
            oprint(msg, dt=(time() - t0), cond=(VERB >= 1))

    del strand_list, bpos_list, epos_list, ovl_idx_cc

# Parallelize the merging of consensus windows if several cores
if (opts['N_PROC'] > 1) and opts['DO_SPOA']:
# if False:
    partial_merge = partial(merge_windows_in_cc, opts=opts)
    pool = Pool(processes=opts['N_PROC'])
    consensi_in_cc = pool.map(partial_merge, range(len(ccs_list)))
    pool.close()
    pool.join()

    for (cc_idx, cons_in_cc) in enumerate(consensi_in_cc):
        print(">contig_%d\n%s" % (cc_idx, cons_in_cc), file=sys.stdout)




oprint("Finished.\nRough layout computed in %4.3f.\n Fine-grained layout computed in %4.3f." % (
        t_rough_layout, t_total_finegrained),
        dt=(time() - t0), cond=(VERB >= 1))

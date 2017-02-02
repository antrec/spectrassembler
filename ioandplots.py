#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""
Functions for Input/Output, including plotting when required.

@author: Antoine Recanati
"""

from __future__ import print_function
import os
import subprocess
import sys
import csv
import matplotlib

matplotlib.use('PS')
import matplotlib.pyplot as plt
import numpy as np


def oprint(msg, dt=None, cond=True):
    """ Prints msg and time elapsed if condition is satisfied.

    Parameters
    ----------
    msg : string (message to print)
    dt : float (time elapsed)
    cond : bool (condition for message to be printed, e.g. verbosity >= 2)

    """
    if not cond:
        return
    elif not dt:
        print(msg, file=sys.stderr)
    else:
        print("%5.3f s: %s" % (dt, msg), file=sys.stderr)


def make_dir(path):
    """ Makes directory given by path if it does not exist yet. """
    if not os.path.exists(path):
        os.mkdir(path)


def fill_args_opts(args):
    """ Checks options from argument parser args in a dictionary and returns them in dictionary opts. """
    # Assign variables and constants and store them in opts dictionary
    opts = {}
    opts['ROOT_DIR'] = args.root
    opts['W_LEN'] = args.w_len
    opts['W_OVL_LEN'] = args.w_ovl_len
    opts['LEN_THR'] = args.len_thr
    opts['VERB'] = VERB = args.verbosity
    opts['MIN_CC_LEN'] = args.min_cc_len
    opts['N_PROC'] = args.nproc
    opts['TRIM_MARGIN'] = args.trim_margin
    opts['MERGE_MARGIN'] = args.margin

    # Check if SPOA found
    if os.path.exists(args.spoapath):
        opts['SPOA_PATH'] = os.path.abspath(args.spoapath)
        opts['DO_SPOA'] = True
    else:
        # In case spoa is in the $PATH but not explicitely given
        try:
            SPOA_PATH = subprocess.check_output("which %s; exit 0" % args.spoapath,
                                                stderr=subprocess.STDOUT, shell=True)
            opts['SPOA_PATH'] = SPOA_PATH.split('\n')[0]
            opts['DO_SPOA'] = True
        # Otherwise do not perform consensus
        except:
            opts['DO_SPOA'] = False
            msg = "spoa executable not found. Provide it with option"\
                  "--spoapath if you wish to compute consensus sequences"
            oprint(msg)

    DO_PLOT_POS_V_REF = False
    if args.ref_pos_csvf is not None:
        DO_PLOT_POS_V_REF = True
        opts['REF_POS_CSVF'] = args.ref_pos_csvf
    opts['DO_PLOT_POS_V_REF'] = DO_PLOT_POS_V_REF

    # Check reads format
    READS_FN = args.READS_FN
    suffix = READS_FN.split('/')[-1].split('.')[-1]
    suffix = suffix[-1]  # only last letter to handle .fasta and .fa the same way
    if (suffix.lower() == 'a'):
        opts['READS_FMT'] = 'fasta'
    elif (suffix.lower() == 'q'):
        opts['READS_FMT'] = 'fastq'
    else:
        msg = "Input file {} has no standard suffix. Please provide a file that ends in "\
              ".*a or .*q (e.g. .fasta, .fa, .fastq or .fq)".format(READS_FN)
        raise StandardError(msg)

    return opts


def plot_cc_pos_v_ref(ref_pos_csvf, cc, bpos_cc, figpath):
    """ Plots position of reads found by algorithm vs reference
    if reference position of the reads is already computed from a reference genome.

    Parameters
    ----------
    all_reads_pos : list - position for all reads in the dataset
    cc : list - reads indices in a given connected component
    cc : list - reads indices in a given connected component
    bpos_cc : list - leftmost base coordinate for reads in cc
    figpath : str - path to file to save figure.

    """
    # Restrict reference position to cc and convert to numpy array
    with open(ref_pos_csvf, "rb") as posf:
        reader = csv.reader(posf)
        all_reads_pos = [int(el) for el in list(reader)[0]]

    ref_pos_cc = np.array(all_reads_pos)[cc]
    bpos_cc = np.array(bpos_cc)
    # Remove reads unmapped by BWA
    ok_idx = np.argwhere(ref_pos_cc < 1e7)
    ref_pos_cc = ref_pos_cc[ok_idx]
    bpos_cc = bpos_cc[ok_idx]

    #      Remove offset
    #    pos_cc -= pos_cc.min()
    #    bpos_cc -= bpos_cc.min()
    plt.scatter(ref_pos_cc, bpos_cc, s=0.1)
    plt.xlabel("true position (from BWA)", fontsize=16)
    plt.ylabel("read position found by algorithm", fontsize=16)
    plt.savefig(figpath)
    plt.close()


def plot_cc_pos_found(bpos_cc, figpath):
    """ For debugging. Plot the fine-grained layout found (positions of reads)
    vs the coarse-grained (ordering of the reads). If it does not look
    *roughly* like a line, there may be a problem.

    Parameters
    ----------
    bpos_cc : list (leftmost base coordinate for reads in current connected component)
    figpath : str (path to file to save figure)

    """
    plt.plot(bpos_cc)
    plt.xlabel("read number (found by spectral ordering)", fontsize=16)
    plt.ylabel("read position found by fine-grained computation", fontsize=16)
    plt.savefig(figpath)
    plt.close()


def write_layout_to_file(layout_fn, strands, bpos, epos, cc, read_nb2id):
    """ Writes layout of reads to a file.

    Writes to file lines containing the following:
    read number (in the order of the fasta file), read name, leftmost base
    coordinate, strand.
    (!) Convention : if strand == '-', the leftmost coordinate corresponds to
    the end of the read (i.e. to the beginning of the reverse component of the
    read).

    Parameters
    ----------
    layout_fn : str (path to file to save layout)
    strands : list (strand for each read)
    bpos : list (leftmost base coordinate for reads current connected component)
    cc : list (read number (from original .fasta/q file) for reads in the current connected component)
    read_nb2id : dict (keys : read number, values : read id)

    """

    idx_sort = bpos.argsort()
    # Write to file
    fh = open(layout_fn, 'wb')
    for idx in idx_sort:
        read_nb = cc[idx]
        read_id = read_nb2id[read_nb]
        lmbc = bpos[idx]
        rmbc = epos[idx]
        s = '+' if strands[idx] else '-'
        fh.write("%d\t%s\t%d\t%d\t%s\n" % (read_nb, read_id, lmbc, rmbc, s))

    fh.close()
    #
    # # Sort reads in connected component by exact position
    # idx_sort = bpos.argsort()
    # bpos = bpos[idx_sort]
    # strands = strands[idx_sort]
    #
    # # Write to file
    # fh = open(layout_fn, 'wb')
    # for (idx, read_nb) in enumerate(cc):
    #     read_id = read_nb2id[read_nb]
    #     lmbc = bpos[idx]
    #     s = '+' if strands[idx] == 1 else '-'
    #     fh.write("%d\t%s\t%d\t%s\n" % (read_nb, read_id, lmbc, s))
    #
    # fh.close()

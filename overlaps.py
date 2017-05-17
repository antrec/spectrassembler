#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""
Tools related to handling overlaps.

Currently implemented to be used with a minimap (https://github.com/lh3/minimap) file.

Includes overlap class, functions that create a sparse matrix from the overlaps,
and a function that computes the absolute positions of the reads from the overlaps in a contig
if the ordering of the reads is given.

@author: Antoine Recanati
"""

import numpy as np
from scipy.sparse import find


class MiniOvl:
    """ Overlap between two reads, named 1 and 2, from line of minimap file.

    Such a line contains :
    query name, length, 0-based start, end, strand,
    target name, length, start, end, the number of matching bases.

    Parameters
    ----------
    mini_line : str (line from minimap file)

    Attributes
    ----------
    id1 : str (read id of read 1)
    id2 : str (read id of read 2)
    len1 : int (length of read 1)
    len2 : int (length of read 2)
    b1 : int (basepair number of the beginning of the overlap on read 1)
    e1 : int (basepair number of the end of the overlap on read 1)
    b2 : int (basepair number of the beginning of the overlap on read 2)
    e2 : int (basepair number of the end of the overlap on read 2)
    strand : char ('+' if the two reads are on same strand and '-' otherwise)
    n_match : int (number of matching bases (see minimap [https://github.com/lh3/minimap] documentation))

    """
    def __init__(self, mini_line):
        fields = mini_line.split()
        self.id1 = fields[0]
        self.len1 = int(fields[1])
        self.b1 = int(fields[2])
        self.e1 = int(fields[3])
        self.strand = fields[4]
        self.id2 = fields[5]
        self.len2 = int(fields[6])
        self.b2 = int(fields[7])
        self.e2 = int(fields[8])
        self.n_match = int(fields[9])
        # self.n_coll = int(fields[10])
        # self.n_frac_match = int(fields[11])

    def switch_ids(self, id1, id2):
        """ Switch reads in the overlap object (read 1 becomes 2 and 2 becomes 1). """
        if (self.id1 == id2) and (self.id2 == id1):
            self.id1, self.id2 = self.id2, self.id1
            self.len1, self.len2 = self.len2, self.len1
            self.b1, self.b2 = self.b2, self.b1
            self.e1, self.e2 = self.e2, self.e1
        else:
            assert self.id1 == id1 and self.id2 == id2, u"id1 : {}, id2 : {} \n self.id1 : {}, self.id2 : {}".format(
                id1, id2, self.id1, self.id2)

    def compute_abs_pos(self, b_ref, s_ref):
        """ Compute absolute position and strand of read 2 from overlap information (self) and absolute
        position and strand of read 1 (b_ref and s_ref).

        Parameters
        ----------
        b_ref : int (absolute position (leftmost base coordinate) of read 1
        s_ref : int (+1 or -1. Absolute strand of read 1)

        Returns
        ----------
        b : int (absolute position (leftmost base coordinate) of read 2)
        s : int (+1 or -1. Absolute strand of read 2)

        """
        # Compute strand of next read
        s = s_ref if self.strand == '+' else not(s_ref)

        # Compute leftmost position (depending of strands of reference and next read)
        if (s_ref and s):
            b = b_ref + self.b1 - self.b2
        elif (s_ref and not(s)):
            b = b_ref + self.b1 - (self.len2 - self.e2)
        elif (not(s_ref) and s):
            b = b_ref + (self.len1 - self.e1) - self.b2
        elif (not(s_ref) and not(s)):
            b = b_ref + (self.len1 - self.e1) - (self.len2 - self.e2)

        return (b, s)


def compute_overlaps(mini_fn, record_list):
    """ Compute list of overlaps from minimap output file and list of reads.

    Parameters
    ----------
    mini_fn : str (path to minimap file)
    record_list : list (list of reads in Bio.SeqIO.records format)

    Returns
    ----------
    read_nb2id : dict (keys : read number, values : read id)
    ovl_list : list (of overlaps as MiniOvl objects)
    i_list : list (of read indices (int) i to build sparse coo_matrix such that A[i,j] ~ overlap between reads i and j)
    j_list : list (of read indices (int) j to build sparse coo_matrix such that A[i,j] ~ overlap between reads i and j)
    k_list : list (of indices (int) k such that ovl_list[k] is the overlap between i_list[k] and j_list[k])
    n_match_list : list (of number of matches (int) such that A[i,j] = number of matches between i and j)
    ovl_len_list : list (of length of overlap between i and j)
    n_reads : int (number of reads)

    """

    # Construct {read name : read number} dictionary
    read_nb_dic = {}
    cpt = 0
    for record in record_list:
        if read_nb_dic.has_key(record.id):
            msg = "Same id {} for reads {} and {} ! " \
                  "Run [https://github.com/antrec/spectrassembler/]check_reads.py "\
                  "on your data first.".format(record.id, read_nb_dic[record.id], cpt)
            raise StandardError(msg)
        read_nb_dic[record.id] = cpt
        cpt += 1
    n_reads = cpt

    idx = 0
    h_list = []
    k_list = []
    ovl_list = []
    n_match_list = []
    ovl_len_list = []
    fh = open(mini_fn, 'rb')

    for line in fh:

        ovl = MiniOvl(line)
        i_idx = read_nb_dic[ovl.id1]
        j_idx = read_nb_dic[ovl.id2]

        # Discard self matches
        if i_idx == j_idx:
            continue

        # Keep 1D indexing : h = n*i + j
        h_idx = n_reads*i_idx + j_idx

        # Check if another overlap between i and j already exists
        duplicate_cond = (h_idx in h_list[-300:])
        if duplicate_cond:
            dupl_idx = h_list[-300:].index(h_idx) + len(h_list) - min(300, len(h_list))
            dupl_ovl = ovl_list[dupl_idx]

            # Drop the overlap if the preexisting one is more significant
            if dupl_ovl.n_match > ovl.n_match:
                continue

            # Replace the preexisting overlap by the new one otherwise
            else:
                n_match_list[dupl_idx] = dupl_ovl.n_match
                ovl_len = (abs(dupl_ovl.e1 - dupl_ovl.b1) \
                    + abs(dupl_ovl.e2 - dupl_ovl.b2))/2
                ovl_len_list[dupl_idx] = ovl_len
                continue

        # Add the overlap if there was no other overlap between i and j
        ovl_list.append(ovl)
        h_list.append(h_idx)
        k_list.append(idx)
        idx += 1
        n_match_list.append(ovl.n_match)
        ovl_len = (abs(ovl.e1 - ovl.b1) + abs(ovl.e2 - ovl.b2))/2
        ovl_len_list.append(ovl_len)

    fh.close()
    # Convert to numpy arrays
    h_list = np.array(h_list)
    n_match_list = np.array(n_match_list)
    ovl_len_list = np.array(ovl_len_list)
    k_list = np.array(k_list)

    # Recover i_list and j_list from h_list indexing (h = n_reads*i + j)
    i_list = h_list//n_reads
    j_list = h_list - n_reads*i_list

    # fh.close()
    read_nb2id = {v : k for (k, v) in read_nb_dic.items()}

    return read_nb2id, ovl_list, i_list, j_list, k_list, n_match_list, ovl_len_list, n_reads


def compute_positions(cc, read_nb2id, ovl_list, ovl_idx_cc):
    """ Compute list of overlaps from minimap output file and list of reads.

    Parameters
    ----------
    cc : list (index of the reads in the cc_idx-th connected component)
    read_nb2id : dict (keys : read number, values : read id)
    ovl_list : list (of overlaps as MiniOvl objects)
    ovl_idx_cc : np.array [n_reads, n_reads] (s.t. ovl_list[ovl_idx_cc[i,j]] = overlap between reads cc[i] and cc[j])

    Returns
    ----------
    strands : numpy.ndarray (of bool, absolute strand for each read in cc
    [+ : True; - : False])
    bpos : numpy.ndarray (of int, absolute leftmost basepair coordinate for each read in cc)
    epos : numpy.ndarray (of int, absolute rightmost basepair coordinate for each read in cc)

    """

    # Count number of neighbors in similarity graph to start with most central
    ovl_ones = ovl_idx_cc.copy()
    ovl_ones.data.fill(1)
    n_nbghr = ovl_ones.sum(axis=0)
    i_start = np.argmax(n_nbghr)

    # Initialization
    cc_len = len(cc)
    strands = np.zeros(cc_len, dtype=bool) # output array containing strand of reads
    bpos = np.zeros(cc_len, dtype=int) # output array containing leftmost coordinate
    #of reads whatever the strand is, i.e. if strand = -1 the leftmost
    #coordinate corresponds to the end of the read
    epos = np.zeros(cc_len, dtype=int) # end of reads (bpos + read_lengths)
    n_nn = 40 # number of neighbors of a read to get the overlaps with
    strands[i_start] = True
    bpos[i_start] = 0
    (_, js, _) = find(ovl_idx_cc[i_start,:])
    ovl_start = ovl_list[ovl_idx_cc[i_start, js[0]]]
    i_id = read_nb2id[cc[i_start]]
    j_id = read_nb2id[cc[js[0]]]
    ovl_start.switch_ids(i_id, j_id)
    epos[i_start] = ovl_start.len1 - 1
    i_done_set = {i_start}
    i_undone_set = set(range(cc_len)) - {i_start}

    ldone_prec = len(i_done_set)
    while len(i_undone_set) > 0:
        # Compute position of reads towards the end of the contig
        for i in xrange(i_start, cc_len):
            (_, i_nbrs, _) = find(ovl_idx_cc[i, :])
            i_nbrs = set(i_nbrs)
            near_nb_set = set(range(i-n_nn/2, i+n_nn/2))
            i_nbrs = i_nbrs.intersection(near_nb_set, i_done_set)
            n_nbrs = len(i_nbrs)
            if (n_nbrs == 0):
                continue
            s_arr = np.zeros(n_nbrs, dtype=bool)
            b_arr = np.zeros(n_nbrs)

            i_id = read_nb2id[cc[i]]
            for (k, j) in enumerate(i_nbrs):
                j_id = read_nb2id[cc[j]]
                ovl_idx = ovl_idx_cc[i, j]
                ovl = ovl_list[ovl_idx]
                ovl.switch_ids(j_id, i_id)
                (b, s) = ovl.compute_abs_pos(bpos[j], strands[j])
                s_arr[k] = s
                b_arr[k] = b

            # Make average
            s = (sum(s_arr) >= 0.5*len(s_arr))
            j_ok = [k for k in range(n_nbrs) if s_arr[k] == s] # read index j that yields the average strand
            b = b_arr[j_ok].mean()
            e = b + ovl.len2 - 1

            strands[i] = s
            bpos[i] = b
            epos[i] = e

            i_done_set = i_done_set | {i}
            i_undone_set = i_undone_set - {i}

        # Compute position of reads towards the beginning of the contig
        for i in xrange(i_start-1, -1, -1):
            (_, i_nbrs, _) = find(ovl_idx_cc[i, :])
            i_nbrs = set(i_nbrs)
            near_nb_set = set(range(i-n_nn/2, i+n_nn/2))
            i_nbrs = i_nbrs.intersection(near_nb_set, i_done_set)
            n_nbrs = len(i_nbrs)
            if (n_nbrs == 0):
                continue
            s_arr = np.zeros(n_nbrs, dtype=bool)
            b_arr = np.zeros(n_nbrs)

            i_id = read_nb2id[cc[i]]
            for (k, j) in enumerate(i_nbrs):
                j_id = read_nb2id[cc[j]]
                ovl_idx = ovl_idx_cc[i, j]
                ovl = ovl_list[ovl_idx]
                ovl.switch_ids(j_id, i_id)
                (b, s) = ovl.compute_abs_pos(bpos[j], strands[j])
                s_arr[k] = s
                b_arr[k] = b

            # Make average
            s = (sum(s_arr) >= 0.5*len(s_arr))
            j_ok = [k for k in range(n_nbrs) if s_arr[k] == s] # read index j that yields the consensus strand
            b = b_arr[j_ok].mean()
            e = b + ovl.len2 - 1

            strands[i] = s
            bpos[i] = b
            epos[i] = e

            i_done_set = i_done_set | {i}
            i_undone_set = i_undone_set - {i}

        # If some reads remain unpositioned in two successive rows,
        # increase the number of neighbors allowed to compute their position
        ldone = len(i_done_set)
        if ldone == ldone_prec:
            n_nn = int(1.5*n_nn)
        ldone_prec = ldone

    return strands, bpos, epos

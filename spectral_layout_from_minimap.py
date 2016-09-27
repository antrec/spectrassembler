# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 10:31:02 2016

@author: Antoine Recanati
"""
from time import time
import os, shutil
import argparse
import numpy as np
from Bio import SeqIO
from scipy.sparse import coo_matrix, eye, find
from scipy.sparse.linalg import eigsh#, ArpackNoConvergence
from scipy.sparse.csgraph import laplacian, connected_components

###############################################################################
############### Functions to handle overlaps from minimap file ################
###############################################################################

class MiniOvl:
    """ Overlap from line of minimap file. By order of appearance :
    query name, length, 0-based start, end, strand,
    target name, length, start, end, the number of matching bases,
    the number of co-linear minimizers in the match
    and the fraction of matching bases. """
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
        self.n_coll = int(fields[10])
        self.n_frac_match = int(fields[11])

def switch_ovl(id1, id2, ovl):
    """ function to make the rest of the code simpler.
    Switch the overlap object ovl between reads id2 and id1
    so as to have ovl.id1 == id1 and ovl.id2 == id2
    If it is already the case, just return the ovl.
    If ovl is not between id1 and id2, raises an error.
    Input : id1, id2, reads identifiers
    ovl, overlap between id1 and id2
    Output : ovl, overlap between id2 and id1 """
    if (ovl.id1 == id1) and (ovl.id2 == id2):
        return ovl

    elif (ovl.id1 == id2) and (ovl.id2 == id1):
        len1 = ovl.len2
        b1 = ovl.b2
        e1 = ovl.e2
        strand = ovl.strand
        len2 = ovl.len1
        b2 = ovl.b1
        e2 = ovl.e1
        n_match = ovl.n_match
        n_coll = ovl.n_coll
        n_frac_match = ovl.n_frac_match
        init_line = "%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d" % \
            (id1, len1, b1, e1, strand, id2, len2, b2, e2, n_match, n_coll, n_frac_match)
        ovlr = MiniOvl(init_line)

        return ovlr
    else:
        raise ValueError("ovl must be overlap between id1 and id2")

def compute_overlaps(mini_fn, record_list):
    """ Build {read name : read number} dictionary,
    list of overlaps (MiniOvl class) ovl_list,
    and overlap index array : ovl_idx_arr[i,j] = k if
    the overlap between reads i and j is ovl_list[k]. """

    # Construct {read name : read number} dictionary
    read_nb_dic = {}
    cpt = 0
    for record in record_list:
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
            dupl_idx = h_list[-300:].index(h_idx) + len(h_list) - 300
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

    fh.close()

    return read_nb_dic, ovl_list, i_list, j_list, k_list, n_match_list, ovl_len_list, n_reads

###############################################################################
############ Functions for pre-processing on the similarity matrix #############
###############################################################################

def sym_max(X):
    """ Function to get maximum between two sparse matrices
    We symmetrize our matrix X with X_sym = max(X, X.T) rather than
    X_sym = X + X.T to avoid adding values in case there are duplicates
    in the overlap file. """
    dif_mat = X - X.T
    dif_mat.data = np.where(dif_mat.data < 0, 1, 0)
    return X - X.multiply(dif_mat) + X.T.multiply(dif_mat)

def remove_bridge_reads(a_mat):
    """ detect "connecting reads" (or bridge reads). i is a connecting read
    if the connected component (of the graph induced by matrix a) that contains i
    is no longer connected if the node i is removed from it. """
    con_reads_list = []
    a_mat = a_mat.copy().tocsr()
    for i in xrange(a_mat.shape[0]):
        (_, J, _) = find(a_mat[i, :])
        if len(J) == 0:
            continue
        Jl = list(set(J))
        a_r = a_mat[Jl, :].tocsc()
        a_r = a_r[:, Jl]
        (n_c, lbl) = connected_components(a_r, directed=False, return_labels=True)
        if n_c >= 2:
            l1 = sum(lbl == 0)
            l2 = sum(lbl == 1)
            if min(l1, l2) > 2:
                con_reads_list.append(i)
    a_clr = a_mat.copy().tolil()
    a_clr[con_reads_list, :] = 0
    a_clr[:, con_reads_list] = 0

    return a_clr

###############################################################################
###### Spectral ordering related functions (gets coarse-grained layout) #######
###############################################################################

def get_fiedler(a):
    """
    Compute the 2nd lowest eigenvalue of the laplacian of a
    and associated eigenvector (Fiedler vector).
    Input : similarity matrix a (must be connected and symmetric.)
    Output : fidval, 2nd smallest eigenvalue
    fidvec : associated eigenvector
    """
    # Construct Laplacian
    lap = laplacian(a, normed=False, return_diag=False) + 1.e-9*eye(a.shape[0])

    # Construct M = lambda_max(lap)*I - lap, whose "second largest"
    # eigenvector is the "second smallest" eigenvector of lap (Fiedler vector)
    (evals_max, _) =  eigsh(lap, 1, which='LA', tol=1e-15, ncv=5)
    maxval = float(evals_max)
    m_lap = maxval*eye(lap.shape[0]) - lap
    if verb >=2:
        print("Computing the Fiedler vector of matrix of size %d ..."\
        % (a.shape[0]))
    evec0 = np.ones((1,lap.shape[0]))
    evals_small, evecs_small = eigsh(m_lap, 2, which='LA', v0=evec0, tol=1e-15, ncv=10)
    if verb >=2:
        print("... Done.")
    fidvec = -evecs_small[:, 0]
    fidval = float(maxval - evals_small[0])

    return fidval, fidvec

def spectral_ordering(a, min_cc_len):
    """ Reorder each connected component (larger than min_cc_len)
    of similarity matrix a with spectral ordering algorithm.
    Input : similarity matrix a,
    min_cc_len (int), minimum length for a connected component to reorder it
    Output : ccs_list, list of reordered connected components.
    ccs_list[k] contains the sequence of reads indices in the connected
    component k reordered by the spectral algorithm.
    n_comp : number of connected components (including those < min_cc_len). """
    ccs_ord = {}
    # Split the similarity graph into connected components \
    #(the spectral algorithms needs a connected similarity graph)
    (n_comp, labels) = connected_components(a, directed=False,
                                            return_labels=True)
    if verb >=2:
        print ("There are %d connected components" % (n_comp))
    for ccn in xrange(n_comp):
        cc = np.argwhere(labels == ccn)[:, 0]
        if len(cc) > min_cc_len:
            # Restrict matrix to the current connected component
            acc = a.tocsr()
            acc = acc[cc, :]
            acc = acc[:, cc]
            # Compute Fiedler vector
            (fidval, fidvec) = get_fiedler(acc)
            # Extract permutation solution and save reordered connected component
            permu = np.argsort(fidvec)
            cc_ord = [cc[idx] for idx in permu]
            ccs_ord[ccn] = cc_ord
            # Check that reordered matrix looks roughly OK
            aloc = a.copy().tocsr()
            aloc = aloc[cc_ord, :]
            aloc = aloc[:, cc_ord]
            (ii, jj, _) = find(aloc)
            bw = max(abs(ii-jj))
            if bw >= 80 and verb >=1:
                print(" WARNING : bandwidth %d >= 80 in reorderd cc %d of size \
                %d. The order found may be wrong" % (bw, ccn, len(cc)))

    # Sort by length of connected component
    ccs_list = ccs_ord.values()
    ccs_len = np.array([len(cc) for cc in ccs_list])
    sort_len = np.argsort(-ccs_len)
    ccs_list = [ccs_list[srtidx] for srtidx in sort_len]

    return ccs_list, n_comp

###############################################################################
############## Functions to compute exact positions of the reads ##############
###############################################################################

def compute_abs_pos(b_ref, s_ref, ovl):
    """ compute absolute position (b : leftmost coordinate of read,
    s = strand compared to first assigned read) from overlap (ovl)
    between already positioned read (begin b_ref strand s_ref) """
    #ovl.id1 : already positioned read; ovl.id2 : read we wish to position
    # Compute strand of next read
    s = s_ref if ovl.strand == '+' else -s_ref

    # Compute leftmost position (depending of strands of reference and next read)
    if ((s_ref == 1) and (s == 1)):
        b = b_ref + ovl.b1 - ovl.b2

    elif ((s_ref == 1) and (s == -1)):
        b = b_ref + ovl.b1 - (ovl.len2 - ovl.e2)

    elif ((s_ref == -1) and (s == 1)):
        b = b_ref + (ovl.len1 - ovl.e1) - ovl.b2

    elif ((s_ref == -1) and (s == -1)):
        b = b_ref + (ovl.len1 - ovl.e1) - (ovl.len2 - ovl.e2)

    return b, s

def compute_positions(cc, read_nb2id, ovl_list, ovl_idx_arr):
    """ Compute strand and leftmost position of reads
    in connected component cc from the overlaps (in ovl_list).
    strands array contains the strand of each read (the first read is assigned
    strand '+'), and bpos array contains the absolute position of the leftmost
    base of the reads (regardless of their orientation).
    epos contains the absolute position of the rightmost coordinates.
    """
    # Restrict overlap index array to reads in the connected component (contig)
    ovl_idx_cc = ovl_idx_arr.tocsc()[:, cc]
    ovl_idx_cc = ovl_idx_cc.tocsr()[cc, :]

    # symmetrize if the overlapper does not count overlap twice for (i,j) and (j,i)
    ovl_idx_cc = sym_max(ovl_idx_cc)

    # Count number of neighbors in similarity graph to start with most central
    ovl_ones = ovl_idx_cc.copy()
    ovl_ones.data.fill(1)
    n_nbghr = ovl_ones.sum(axis=0)
    i_start = np.argmax(n_nbghr)

    # Initialization
    cc_len = len(cc)
    strands = np.empty(cc_len, dtype=int) # output array containing strand of reads
    bpos = np.empty(cc_len, dtype=int) # output array containing leftmost coordinate
    #of reads whatever the strand is, i.e. if strand = -1 the leftmost
    #coordinate corresponds to the end of the read
    epos = np.empty(cc_len, dtype=int) # end of reads (bpos + read_lengths)
    n_nn = 40 # number of neighbors of a read to get the overlaps with
    strands[i_start] = 1
    bpos[i_start] = 0
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
            s_arr = np.empty(n_nbrs)
            b_arr = np.empty(n_nbrs)

            i_id = read_nb2id[cc[i]]
            for (k, j) in enumerate(i_nbrs):
                j_id = read_nb2id[cc[j]]
                ovl_idx = ovl_idx_cc[i, j]
                ovl = ovl_list[ovl_idx]
                ovl = switch_ovl(j_id, i_id, ovl)
                (b, s) = compute_abs_pos(bpos[j], strands[j], ovl)
                s_arr[k] = s
                b_arr[k] = b

            # Make average
            s = np.sign(sum(s_arr)) #strand
            j_ok = [k for k in range(n_nbrs) if s_arr[k] == s] # read index j that yields the average strand
            b = np.mean(b_arr[j_ok])
            e = b + ovl.len2

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
            s_arr = np.empty(n_nbrs)
            b_arr = np.empty(n_nbrs)

            i_id = read_nb2id[cc[i]]
            for (k, j) in enumerate(i_nbrs):
                j_id = read_nb2id[cc[j]]
                ovl_idx = ovl_idx_cc[i, j]
                ovl = ovl_list[ovl_idx]
                ovl = switch_ovl(j_id, i_id, ovl)
                (b, s) = compute_abs_pos(bpos[j], strands[j], ovl)
                s_arr[k] = s
                b_arr[k] = b

            # Make average
            s = np.sign(sum(s_arr)) #strand
            j_ok = [k for k in range(n_nbrs) if s_arr[k] == s] # read index j that yields the consensus strand
            b = np.mean(b_arr[j_ok])
            e = b + ovl.len2

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

###############################################################################
############## Functions to write layout and input files for POA ##############
###############################################################################

def write_layout_to_file(layout_fn, strands, bpos, cc, read_nb2id):
    """ writes layout of reads to a file named layout_fn. Each line contains :
    read number (in the order of the fasta file), read name, leftmost base
    coordinate, strand.
    (!) Convention : if strand == '-', the leftmost coordinate corresponds to
    the end of the read (i.e. to the beginning of the reverse component of the
    read). """

    # Sort reads in connected component by exact position
    idx_sort = np.argsort(bpos)
    bpos = bpos[idx_sort]
    strands = strands[idx_sort]

    # Write to file
    fh = open(layout_fn, 'wb')
    for (idx, read_nb) in enumerate(cc):
        read_id = read_nb2id[read_nb]
        lmbc = bpos[idx]
        s = '+' if strands[idx] == 1 else '-'
        fh.write("%d\t%s\t%d\t%s\n" % (read_nb, read_id, lmbc, s))

    fh.close()
    return

def write_poa_inputs(record_list, cc_dir, cc_idx, cc, strand_list, bpos_list,
                     epos_list, read_nb2id, w_len, w_ovl_len):
    """ Writes input files for POA for a given connected component (cc).
    Each file contains the portions of reads sequences that fit in a given
    window.
    record_list : list of reads/records (SeqIO format)
    cc_dir : directory where to write the files for this cc
    cc_idx : which cc
    cc : oredered sequence of reads in the connected component number cc_idx
    strand_list, bpos_list, epos_list : strand, leftmost and rightmost
    coordinates of the reads ordered as in cc
    read_nb2id : dictionary {read index : read identifier}
    w_len : length of the windows
    w_ovl_len : length of the overlap between two consecutive windows"""

    n_reads = len(bpos_list)
    idx_range_list = range(n_reads)
    if not(bpos_list.min() == 0):
        min_b = bpos_list.min()
        bpos_list -= min_b
        epos_list -= min_b

    cons_len = epos_list.max()
    n_windows = cons_len // (w_len - w_ovl_len)
    n_windows = int(n_windows +1)

    for w_idx in xrange(n_windows):
        w_b = (w_len - w_ovl_len)*w_idx
        w_e = w_b + w_len
        reads_in_w = [idx for idx in idx_range_list if (bpos_list[idx] <= w_e and epos_list[idx] >= w_b)]
        in_fn = cc_dir + "/poa_in_cc_%d_win_%d.fasta" % (cc_idx, w_idx)
        in_fh = open(in_fn, 'wb')
        for idx in reads_in_w:
            read_idx = cc[idx]
            read_id = read_nb2id[read_idx]
            record = record_list[read_idx]
            read_seq = record.seq
            # Reverse read if strand is '-'
            strd = strand_list[idx]
            if strd == -1:
                read_seq = read_seq.reverse_complement()

            # Trim read to the part contained in the window
            read_len = len(read_seq)
            bb = int(max(0, w_b - bpos_list[idx]))
            ee = int(read_len - max(0, epos_list[idx] - w_e))
            read_seq = read_seq[bb:ee]

            # Write to poa_in file
            in_fh.write(">%s\n%s\n" % (read_id, read_seq))

        in_fh.close()

    return

###############################################################################
############### Functions for plotting (if reference available ################
###############################################################################

def plot_cc_pos_v_ref(all_reads_pos, cc, bpos_cc, figpath):
    """ If reference genome available and position of the reads computed.
    all_reads_pos : list of position for all reads in the dataset
    cc : list of reads indices in a given connected component
    bpos_cc : leftmost base coordinate for reads in cc
    figpath : path to file to save figure. """
    # Restrict reference position to cc and convert to numpy array
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
    return

def plot_cc_pos_found(bpos_cc, figpath):
    """ Sanity check : plot the fine-grained layout found (positions of reads)
    vs the coarse-grained (ordering of the reads). If it does not look
    *roughly* like a line, there may be a problem.
    bpos_cc : leftmost base coordinate for reads in current connected component
    figpath : path to file to save figure. """
    plt.plot(bpos_cc)
    plt.xlabel("read number (found by spectral ordering)", fontsize=16)
    plt.ylabel("read position found by fine-grained computation", fontsize=16)
    plt.savefig(figpath)
    plt.close()
    return

###############################################################################
###############################################################################
###############################################################################
#################################### Main #####################################
###############################################################################
###############################################################################
###############################################################################

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--root", type=str, default="./",
                    help="directory where to store layout files.")
parser.add_argument("-f", "--fastafn", type=str, required=True,
                    help="reads file path (fasta format))")
parser.add_argument("-m", "--minimapfn", type=str, required=True,
                    help="overlap file path (from minimap in PAF format).")
parser.add_argument("-w", "--write_poa_files", action="store_true",
                    help="Whether to write POA input files for \
consensus generation or not.")
parser.add_argument("--w_len", type=int, default=2500,
                    help="length of consensus windows for POA.")
parser.add_argument("--w_ovl_len", type=int, default=1250,
                    help="overlap length between two successive \
consensus windows.")
parser.add_argument("--sim_thr", type=int, default=850,
                    help="threshold on overlap score \
(similarity matrix preprocessing.)")
parser.add_argument("--len_thr", type=int, default=3500,
                    help="threshold on length of overlaps \
(similarity matrix preprocessing).")
parser.add_argument("-v", "--verbosity", action="count", default=0,
                    help="verbosity level (-v, -vv or none)")
parser.add_argument("--ref_pos_csvf",
                    help="csv file with position of reads \
(in same order as in fastafn) obtained from BWA, \
in order to plot reads position found vs reference.")

args = parser.parse_args()

# Assign variables and constants
root_dir = args.root
w_len = args.w_len
w_ovl_len = args.w_ovl_len
THR = args.sim_thr
LEN_THR = args.len_thr
do_write_poa = args.write_poa_files
global verb
verb = args.verbosity
do_plot_pos = False
#do_plot_pos = True if verb >=2 else False
do_plot_pos_v_ref = False
if not(args.ref_pos_csvf == None):
    do_plot_pos_v_ref = True
    ref_pos_csvf = args.ref_pos_csvf
    import csv
    with open(ref_pos_csvf, "rb") as posf:
        reader=csv.reader(posf)
        ref_pos_list = [ int(el) for el in list(reader)[0] ]

if (do_plot_pos or do_plot_pos_v_ref):
    import matplotlib
    matplotlib.use('PS')
    import matplotlib.pyplot as plt

# Start main program
t0 = time()
t1 = t0
# Load fasta file
if verb >= 2:
    print("Loading reads...")
reads_fh = open(args.fastafn, "rU")
record_list = list(SeqIO.parse(reads_fh, "fasta"))
reads_fh.close()
if verb >=2:
    print("... Done in %3.1fs" % (time() - t1))
    t1 = time()
    print("Compute overlaps from files...")
# Compute overlaps from the files
(read_id2idx, ovl_list, I, J, K, num_match, ovl_len, n_reads) = \
    compute_overlaps(args.minimapfn, record_list)

if verb >= 2:
    print("... Done in %3.1fs" % (time() - t1))

# Revert read_id to read_idx dictionary
read_nb2id = {v : k for (k, v) in read_id2idx.items()}

### Preprocessing and building of the similarity matrix ###
# Threshold based on overlaps value (number of matches) and length
cond1 = (num_match > THR)
cond2 = (ovl_len > LEN_THR)
idxok = np.argwhere(cond1*cond2)[:, 0]
I = I[idxok]
J = J[idxok]
num_match = num_match[idxok]
ovl_len = ovl_len[idxok]
K = K[idxok]

# Construct the matrix
if verb >= 2:
    print("Construct thresholded similarity matrix...")
    t1 = time()
sim_mat = coo_matrix((num_match, (I, J)), shape=(n_reads, n_reads), dtype=int)
if verb >=2:
    print("... Done in %3.1fs" % (time() - t1))
    t1 = time()
    print("Preprocess similarity matrix...")
# Overlap index array : overlap(i,j) = ovl_list[k], with k = ovl_idx_arr[i,j]
ovl_idx_arr = coo_matrix((K, (I, J)), shape=(n_reads, n_reads), dtype=int)

# Symmetrize the matrix when it is not already symmetric
sim_mat = sym_max(sim_mat)
#sim_mat = (sim_mat + sim_mat.T)

# Remove "connecting reads"
sim_mat = remove_bridge_reads(sim_mat)

if verb == 1:
    print ("Similarity matrix built and preprocessed. \
    Total time elapsed : %3.1fs." % (time() - t0))
elif verb >=2:
    print("... Done in %3.1fs" % (time() - t1))
    t1 = time()
if verb >=2:
    print ("Reorder similarity matrix with spectral ordering.")
    t1 = time()
### Compute the layout of the reads ###
# Reorder connected components with spectral ordering
MIN_CC_LEN = 3 # minimum number of reads in a connected component to layout
(ccs_list, n_comp) = spectral_ordering(sim_mat, MIN_CC_LEN)
if verb == 1:
    print ("Similarity matrix reordered. Total time elapsed %3.1fs" \
    % (time() - t0))
elif verb >=2:
    print ("... Done in %3.1fs. Total time elapsed %3.1fs" \
    % (time() - t1, time() - t0))

# Get fine-grained layout with dictionary of overlaps in each connected component
if verb >= 1:
    print("Compute fine-grained layout, print positions of reads to file \
    and print POA input files and graphics if asked for.")
t1 = time()

# If root_dir does not exist, create it
try:
    os.mkdir(root_dir)
    os.chdir(root_dir)
except:
    os.chdir(root_dir)

for (cc_idx, cc) in enumerate(ccs_list):

    # Compute strand and leftmost base coordinate for each read in cc
    (strand_list, bpos_list, epos_list) = compute_positions(cc, read_nb2id, ovl_list, ovl_idx_arr)
    if verb >=2:
        print("Connected component %d/%d" % (cc_idx, len(ccs_list)))
        print "Positions computed in %3.1fs - cc %d of size %d" % (time() - t1, cc_idx +1, len(cc))

    # Write file with layout
    layout_fn = "%s/cc%d.layout" % (root_dir, cc_idx)
    write_layout_to_file(layout_fn, strand_list, bpos_list, cc, read_nb2id)
    if verb >= 2:
        print("Layout written to file.")
        t1 = time()

    if do_write_poa:
        # Write input files for POA to compute consensus
        cc_dir = "%s/cc_%d" % (root_dir, cc_idx)
        try:
            os.mkdir(cc_dir)
        except:
            print "directory for POA input files already exists. Its content will be deleted and rewritten."
            shutil.rmtree(cc_dir)
            os.mkdir(cc_dir)
        write_poa_inputs(record_list, cc_dir, cc_idx, cc, strand_list, bpos_list, epos_list, read_nb2id, w_len, w_ovl_len)
        if verb >= 2:
            print("Written POA input files to %s/cc_%d/." % (root_dir, cc_idx))
            t1 = time()

    if do_plot_pos:
        if verb >= 2:
            print("Edit graphic (position of reads found by algorithm).")
        figpath = root_dir + "/pos_found_cc%d.eps" % (cc_idx)
        plot_cc_pos_found(bpos_list, figpath)

    if do_plot_pos_v_ref:
        if verb >= 2:
            print("Edit graphic (position of reads found by algorithm vs \
            reference (BWA)).")
        figpath = root_dir + "/pos_found_vs_ref_cc%d.eps" % (cc_idx)
        plot_cc_pos_v_ref(ref_pos_list, cc, bpos_list, figpath)

if verb >= 1:
    print("All done, total time elapsed : %3.1fs." % (time() - t0))

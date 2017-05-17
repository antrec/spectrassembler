#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tools related to the spectral ordering algorithm.

Includes functions for pre-processing of a sparse similarity matrix,
and a recursive function that reorders its connected component and keeps on
increasing a threshold on the values of the similarity matrix if the ordering found
does not exhibit a sufficiently small bandwidth.

@author: Antoine Recanati
"""
from time import time
import subprocess
import os
# import sys
from functools import partial
from itertools import repeat
from multiprocessing import Pool
import numpy as np
from scipy.sparse import coo_matrix, eye, find, isspmatrix_csr
from scipy.sparse.linalg import eigsh  # , ArpackNoConvergence
from scipy.sparse.csgraph import laplacian, connected_components
from scipy.stats.mstats import mquantiles
from random import randrange

from ioandplots import oprint


### Functions for pre-processing

def sym_max(X):
    """
    Returns symmetrization of sparse matrix X.

    X_sym = max(X, X.T) rather than X + X.T to avoid adding up values when there are duplicates in the overlap file.
    If X is triangular, max(X, X.T) and X + X.T are equal. """

    dif_mat = X - X.T
    dif_mat.data = np.where(dif_mat.data < 0, 1, 0)
    return X - X.multiply(dif_mat) + X.T.multiply(dif_mat)


def remove_bridge_reads(a_mat):
    """ Remove some edges from the similarity graph.

    When the set of neighbors N(i) of a node i are not connected if that node i is removed from the graph,
    the edges between i and j are cut for all j that is not in the largest connected group among N(i).

    Parameters
    ----------
    a_mat : scipy.sparse matrix (similarity matrix)

    Returns
    ----------
    a_clr : scipy.sparse matrix (similarity matrix pre-preocessed)

    """
    Ikill = []
    Jkill = []
    if not(isspmatrix_csr(a_mat)):
        a_mat = a_mat.tocsr()
    for i in xrange(a_mat.shape[0]):
        (_, J, _) = find(a_mat[i, :])
        if len(J) == 0:
            continue
        Jl = list(set(J))
        a_r = a_mat[Jl, :].tocsc()
        a_r = a_r[:, Jl]
        Jl = np.array(Jl)
        (n_c, lbl) = connected_components(a_r, directed=False, return_labels=True)
        if n_c > 1:
            sizeccs = np.zeros(n_c)
            for ncc in xrange(n_c):
                sizeccs[ncc] = sum(lbl == ncc)
            ccmax = np.argmax(sizeccs)
            away_idx = np.where(lbl != ccmax)[0]
            away_nbrs = list(Jl[away_idx])
            Ikill.extend([i] * len(away_nbrs))
            Jkill.extend(away_nbrs)

    Ikill = np.array(Ikill)
    Jkill = np.array(Jkill)
    Vkill = np.ones(Ikill.size)
    kill_mat = coo_matrix((Vkill, (Ikill, Jkill)), shape=a_mat.shape, dtype=int).tocsr()
    kill_mat = sym_max(kill_mat)
    kill_mat = kill_mat.multiply(a_mat)
    a_clr = a_mat - kill_mat
    if not(isspmatrix_csr(a_clr)):
        a_clr = a_clr.tocsr()

    return a_clr


###############################################################################
###### Spectral ordering related functions (gets coarse-grained layout) #######
###############################################################################

def get_fiedler(A):
    """
    Compute the 2nd lowest eigenvalue of the laplacian of A and associated eigenvector (Fiedler vector).

    Parameters
    ----------
    A : scipy.sparse matrix (similarity matrix)

    Returns
    ----------
    fidval : real (2nd smallest eigenvalue)
    fidvec : array (associated eigenvector)

    """
    # Construct Laplacian
    L_A = laplacian(A, normed=False, return_diag=False) + 1.e-9 * eye(A.shape[0])

    # Construct M = lambda_max(lap)*I - lap, whose "second largest"
    # eigenvector is the "second smallest" eigenvector of lap (Fiedler vector)
    (evals_max, _) = eigsh(L_A, 1, which='LA', tol=1e-15, ncv=20)
    maxval = float(evals_max)
    m_lap = maxval * eye(L_A.shape[0]) - L_A
    evec0 = np.ones((1, L_A.shape[0]))
    evals_small, evecs_small = eigsh(m_lap, 2, which='LA', v0=evec0, tol=1e-20, ncv=20)
    fidvec = -evecs_small[:, 0]
    fidval = float(maxval - evals_small[0])

    return fidval, fidvec

def get_fiedler_julia(mat, julia_path, julia_fiedler_script):

    # Construct Laplacian
    (iis, jjs, vvs) = find(mat)
    n = mat.shape[0]
    # add a random number to the filename to avoid inconsistent writes and reads with multiprocessing
    randn = randrange(1000)

    # write to csv
    itempf = 'mat_coords_iis_%d_%d.csv' % (n, randn)
    iis.tofile(itempf, sep=',')
    jtempf = 'mat_coords_jjs_%d_%d.csv' % (n, randn)
    jjs.tofile(jtempf, sep=',')
    vtempf = 'mat_data.csv_%d_%d' % (n, randn)
    vvs.tofile(vtempf, sep=',')

    outf = 'temp_fidvec_%d_%d.csv' % (n, randn)
    # call julia
    cmd = [julia_path, julia_fiedler_script, itempf, jtempf, vtempf, outf]
    subprocess.call(cmd)

    # remove temporary files
    os.remove(itempf)
    os.remove(jtempf)
    os.remove(vtempf)

    # check output looks OK and return permutation
    if os.path.exists(outf):
        myperm = np.fromfile(outf, dtype=int, sep=',')
        myperm = myperm - 1
        os.remove(outf)
        if (len(myperm) == mat.shape[0]):
            return myperm
        else:
            return np.arange(n)
    # output identity permutation if something went wrong
    else:
        return np.arange(n)

def reorder_submat(A, cc, num_match_l, qtile, ccs_ord, opts):
    """ Reorder matrix A with spectral ordering algorithm.

    Recursive function that reorders each connected component of the input matrix and raises threshold in
    the connected components where the order seems wrong, based on the bandwidth of the reordered matrix
    (this criterium is empirical and specific to genome assembly of genomes with limited number of repeats).

    Parameters
    ----------
    A : scipy.sparse matrix (similarity matrix)
    cc : list (index of the reads in the cc_idx-th connected component)
    num_match_l : list (of number of matches (int) such that A[i,j] = number of matches between i and j) *before*
    preprocessing and not restricted to the reads in cc. It is used to compute the threshold with qtile)
    qtile : real (the values lower than the threhsold thr = quantile(num_match_l, qtile) are removed from A)
    opts : dict (keywords argument containing global parameters and options)

    ccs_ord : list (of lists or reads index sorted by position inside a given connected component)


    Returns
    ----------
    None but ccs_ord is modified "passed by reference"

    """

    VERB = opts['VERB']
    min_cc_len = opts['MIN_CC_LEN']
    JULIA_PATH = opts['JULIA_PATH']
    JULIA_SCRIPT = opts['JULIA_SCRIPT']

    # rep_time_fh = open('%s/time_evs.txt' %(opts['ROOT_DIR']), 'wb')
    # t0 = time()

    if not isspmatrix_csr(A):
        A = A.tocsr()
    (ncs, lbls) = connected_components(A, directed=False, return_labels=True)
    for nc in xrange(ncs):
        cc_sub = np.argwhere(lbls == nc)[:, 0]
        if len(cc_sub) <= min_cc_len:
            continue
        msg = " Running spectral algorithm in connected component of size %d..." % (len(cc_sub))
        oprint(msg, cond=(VERB >= 2))
        # A_sub = A.copy().tocsr()
        # A_sub = A_sub[cc_sub, :]
        # A_sub = A_sub[:, cc_sub]
        A_sub = A[cc_sub, :][:, cc_sub]
        # t1 = time()
        #

        # Use Julia if possible to reorder relatively large matrices
        if JULIA_PATH and (len(cc_sub) > 4000):
            permu = get_fiedler_julia(A_sub, JULIA_PATH, JULIA_SCRIPT)
            # rep_time_fh.write("%d\t%3.6f\t(julia)\n" %(len(cc_sub), time()-t1))

        else:
            (fidval, fidvec) = get_fiedler(A_sub)
            if fidval < 1e-12:
                oprint("\n\nWARNING ! Non connected submatrix of size %d!\n\n" % (len(cc_sub)))
            # rep_time_fh.write("%d\t%3.6f\n" %(len(cc_sub), time()-t1))
            permu = np.argsort(fidvec)
        cc_ord = [cc_sub[idx] for idx in permu]
        # A_ord = A_sub.copy()
        # A_ord = A_ord[permu, :]
        # A_ord = A_ord[:, permu]
        # (ii, jj, _) = find(A_ord)
        (ii, jj, _) = find(A_sub[permu, :][:, permu])
        bw = max(abs(ii - jj))
        if bw >= 80:
            oprint("Bandwidth larger than 80 in reordered matrix. Threshold in submatrix increased before reordering.",
                   cond=(VERB >= 2))
            new_qtile = qtile
            new_qtile += min(0.1, 0.5 * (1. - qtile))
            thr_sub = mquantiles(num_match_l, new_qtile)
            A_sub = remove_bridge_reads(A_sub.multiply(A_sub > thr_sub))
            cc_abs = [cc[idx] for idx in cc_sub]
            reorder_submat(A_sub, cc_abs, num_match_l, new_qtile, ccs_ord, opts)
        else:
            ccs_ord.append([cc[idx] for idx in cc_ord])
            # oprint("Done in %3.3f." %(time() - t1), dt=(time() - t0), cond=(VERB >= 2))
        #
        # oprint("Computed rough layout in %3.3f." %(time() - t0), cond=(VERB >= 2))
        #
    # rep_time_fh.close()
    return

def reorder_mat(A, thr_list, min_cc_len, VERB):

    if not isspmatrix_csr(A):
        A = A.tocsr()
    # Initialization.
    ccs_ord = []
    #Create list of unordered connected components
    todo_ccs = [np.arange(A.shape[0])]
    todo_next = []
    n_loop = 0

    while len(todo_ccs) > 0:
        thr_sub = thr_list[n_loop] # starts at 0.4 for n_loop=0
        # Reorder each of them
        for cc in todo_ccs:
            # if statement
            # in order not to make the preprocessing twice. We could also remove
            # the preprocessing from the pipeline and do it here.
            if n_loop > 0:
                A_sub = A[cc, :][:, cc]
                A_sub = remove_bridge_reads(A_sub.multiply(A_sub > thr_sub))
            else:
                A_sub = A

            # Compute connected components
            (n_cc, labels) = connected_components(A_sub, directed=False, return_labels=True)

            # Reorder each cc with spectral and keep the ordering if it looks OK
            for i_cc in xrange(n_cc):
                cc_sub = np.argwhere(labels == i_cc)[:, 0]
                if len(cc_sub) <= min_cc_len:
                    continue
                msg = " Running spectral algorithm in connected"\
                      "component of size %d..." % (len(cc_sub))
                oprint(msg, cond=(VERB >= 2))
                (_, fidvec) = get_fiedler(A_sub[cc_sub, :][:, cc_sub])
                permu = np.argsort(fidvec)
                (ii, jj, _) = find(A_sub[cc_sub[permu], :][:, cc_sub[permu]])
                bw = max(abs(ii - jj))
                if bw >= 80:
                    oprint("Bandwidth larger than 80 in reordered matrix.",
                           cond=(VERB >= 2))
                    todo_next.append(cc[cc_sub])
                else:
                    ccs_ord.append(cc[cc_sub[permu]])

        todo_ccs = todo_next
        todo_next = []
        n_loop += 1

    return ccs_ord


def reord_submat(in_tuple, A, opts):

    (thr_sub, cc) = in_tuple
    min_len = int(opts['MIN_CC_LEN'])
    verb = int(opts['VERB'])
    JULIA_PATH = opts['JULIA_PATH']
    JULIA_SCRIPT = opts['JULIA_SCRIPT']
    # rep_time_fh = open('%s/time_evs.txt' %(opts['ROOT_DIR']), 'a')

    sub_todo_next = []
    sub_ccs_ord = []

    A_sub = A[cc, :][:, cc]
    A_sub = remove_bridge_reads(A_sub.multiply(A_sub > thr_sub))
    # Compute connected components
    (n_cc, labels) = connected_components(A_sub, directed=False, return_labels=True)

    # Reorder each cc with spectral and keep the ordering if it looks OK
    for i_cc in xrange(n_cc):
        cc_sub = np.argwhere(labels == i_cc)[:, 0]
        if len(cc_sub) <= min_len:
            continue
        msg = " Running spectral algorithm in connected "\
              "component of size %d..." % (len(cc_sub))
        oprint(msg, cond=(verb >= 2))
        # t1 = time()

        if JULIA_PATH and (len(cc_sub) > 3000):
            permu = get_fiedler_julia(A_sub[cc_sub, :][:, cc_sub], JULIA_PATH, JULIA_SCRIPT)
            # rep_time_fh.write("%d\t%3.6f\t(julia)\n" %(len(cc_sub), time()-t1))
        else:
            (_, fidvec) = get_fiedler(A_sub[cc_sub, :][:, cc_sub])
            permu = np.argsort(fidvec)
            # rep_time_fh.write("%d\t%3.6f\n" %(len(cc_sub), time()-t1))

        oprint("Done in %3.6fs" % (time()-t1), cond=(verb>=2))

        (ii, jj, _) = find(A_sub[cc_sub[permu], :][:, cc_sub[permu]])
        bw = max(abs(ii - jj))
        if bw >= 80:
            oprint("Bandwidth larger than 90 in reordered matrix.",
                   cond=(verb >= 2))
            sub_todo_next.append(cc[cc_sub])
        else:
            sub_ccs_ord.append(cc[cc_sub[permu]])

    # rep_time_fh.close()

    return sub_ccs_ord, sub_todo_next

def reorder_mat_par(A, thr_list, opts):

    partial_reorder = partial(reord_submat, A=A, opts=opts)
    N_PROC = int(opts['N_PROC'])//4
    min_cc_len = opts['MIN_CC_LEN']

    if not isspmatrix_csr(A):
        A = A.tocsr()
    # Initialization.
    ccs_ord = []
    #Create list of unordered connected components
    todo_ccs = [np.arange(A.shape[0])]
    todo_next = []
    n_loop = 0

    todo_ccs = []
    (ncs, lbls) = connected_components(A, directed=False, return_labels=True)
    for nc in xrange(ncs):
        cc_sub = np.argwhere(lbls == nc)[:, 0]
        if len(cc_sub) <= min_cc_len:
            continue
        todo_ccs.append(cc_sub)


    while len(todo_ccs) > 0:
        thr_sub = thr_list[n_loop] # starts at 0.4 for n_loop=0
        args_list = zip(repeat(thr_sub), todo_ccs)

        pool = Pool(processes=N_PROC)
        results = pool.map(partial_reorder, args_list)
        pool.close()
        pool.join()

        for tple in results:
            (sub_ccs_ord, sub_todo_next) = tple
            ccs_ord += sub_ccs_ord
            todo_next += sub_todo_next

        todo_ccs = todo_next
        todo_next = []
        n_loop += 1
    return ccs_ord

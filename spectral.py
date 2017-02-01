# -*- coding: utf-8 -*-
"""
Tools related to the spectral ordering algorithm.

Includes functions for pre-processing of a sparse similarity matrix,
and a recursive function that reorders its connected component and keeps on
increasing a threshold on the values of the similarity matrix if the ordering found
does not exhibit a sufficiently small bandwidth.

@author: Antoine Recanati
"""
import warnings
import numpy as np
from scipy.sparse import coo_matrix, eye, find
from scipy.sparse.linalg import eigsh  # , ArpackNoConvergence
from scipy.sparse.csgraph import laplacian, connected_components
from scipy.stats.mstats import mquantiles

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
    kill_mat = coo_matrix((Vkill, (Ikill, Jkill)), shape=a_mat.shape, dtype=int)
    kill_mat = sym_max(kill_mat)
    kill_mat = kill_mat.multiply(a_mat)
    a_clr = a_mat - kill_mat

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
    (evals_max, _) = eigsh(L_A, 1, which='LA', tol=1e-30, ncv=5)
    maxval = float(evals_max)
    m_lap = maxval * eye(L_A.shape[0]) - L_A
    evec0 = np.ones((1, L_A.shape[0]))
    evals_small, evecs_small = eigsh(m_lap, 2, which='LA', v0=evec0, tol=1e-30, ncv=15)
    fidvec = -evecs_small[:, 0]
    fidval = float(maxval - evals_small[0])

    return fidval, fidvec


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

    (ncs, lbls) = connected_components(A, directed=False, return_labels=True)
    for nc in xrange(ncs):
        cc_sub = np.argwhere(lbls == nc)[:, 0]
        if len(cc_sub) <= min_cc_len:
            continue
        msg = " Running spectral algorithm in connected component of size %d..." % (len(cc_sub))
        oprint(msg, cond=(VERB >= 2))
        A_sub = A.copy().tocsr()
        A_sub = A_sub[cc_sub, :]
        A_sub = A_sub[:, cc_sub]
        (fidval, fidvec) = get_fiedler(A_sub)
        if fidval < 1e-12:
            oprint("\n\nWARNING ! Non connected submatrix of size %d!\n\n" % (len(cc_sub)))
        permu = np.argsort(fidvec)
        cc_ord = [cc_sub[idx] for idx in permu]
        A_ord = A_sub.copy()
        A_ord = A_ord[permu, :]
        A_ord = A_ord[:, permu]
        (ii, jj, _) = find(A_ord)
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
            oprint("Done.", cond=(VERB >= 2))
    return

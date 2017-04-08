# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 16:25:53 2015

@author: Antoine Recanati
"""

import argparse
import time
import numpy as np
import pandas as pd
from Bio import SeqIO

import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt


# reverse complement a sequence
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def revcomp(seq):
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement

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


def make_fragment_maps(fragment, RE_SITE, minLen):
    try:
        L =   [len(el) + len(RE_SITE) for el in fragment.split(RE_SITE)[1:-1]]
        Lrev = [len(el) + len(RE_SITE) for el in revcomp(fragment).split(revcomp(RE_SITE))[1:-1]]
        L2 = []; L2rev = []
        tmp = 0
        for i in range(len(L)):
            tmp += L[i]
            if tmp > minLen:
                L2.append(tmp)
                tmp = 0
        tmp = 0
        for i in range(len(Lrev)):
            tmp += Lrev[i]
            if tmp > minLen:
                L2rev.append(tmp)
                tmp = 0
        beginLength = len(fragment.split(RE_SITE)[0])
        endLength = len(fragment.split(RE_SITE)[-1])
        return L2, L2rev, beginLength, endLength
    except:
        print "not enough RE occurences"
        return [],[], 0, 0

def find_best_match_DP(fragMap, refMap, C, sigmaInv, addLength, maxShift):
    """
    Dynamic Programming algorithm as described in the paper
    ''Scaffolding and validation of bacterial genome assemblies
    using optical restriction maps''
    of  Niranjan Nagarajan, Timothy D. Read and Mihai Pop
    """
    frag_map = fragMap[:]
    len_frag = len(frag_map)
    ref_map_loc = refMap[:]
    len_ref = len(ref_map_loc)
    S = np.zeros((len_frag+1, len_ref+1))
    S[:, :] = -float("inf")
    S[0, :] = 0.
    for i in range(1, len_frag+1):
        for j in range(1, len_ref+1):
            iMin = max(0, i-maxShift)
            jMin = max(0, j-maxShift)
            (xx, yy) = np.meshgrid(range(i-iMin-1, -1, -1), range(j-jMin-1, -1, -1))
            (cc, oo) = np.meshgrid(np.cumsum(frag_map[iMin:i][::-1])[::-1], \
            np.cumsum(ref_map_loc[jMin:j][::-1])[::-1])
            (_, nn) = np.meshgrid(np.ones_like(range(i-iMin)), range(j-jMin, 0, -1))
            devMat = -C*abs((xx + yy)).T
            distMat = - sigmaInv*(np.divide(abs(cc - oo), nn)).T
            prevMat = S[iMin:i, jMin:j]
            subMat = (devMat + distMat + prevMat).copy()
            S[i, j] = subMat.max()
    jm = np.argmax(S[-1, :])
    score = S[-1, jm]
    # TEST AU 18/09/2016 : refMap[:jm + 1]
    # endPos = sum(refMap[:jm]) + addLength
    endPos = sum(refMap[:jm+1]) + addLength
    return endPos, score

# Dictionary of Restriction Enzymes
enz_dic = {'AflII' : 'CTTAAG', 'BamHI' : 'GGATCC', 'EcoRI' : 'GAATTC', \
           'EcoRV' : 'GATATC', 'NcoI' : 'CCATGG', 'NotI' : 'GCGGCCGC', \
           'NruI' : 'TCGCGA', 'PvuII' : 'CAGCTG', 'SalI' : 'GTCGAC', \
           'ScaI' : 'AGTACT', 'XbaI' : 'TCTAGA', 'XhoI' : 'CTCGAG'}

parser = argparse.ArgumentParser(description="Simulates optical mapping experiment"
                                             "with contigs from S288C yeast genome"
                                             "and a perfect optical map simulated in silico")
parser.add_argument("sequences",
                    help="sequences output from spectrassembler in FASTA format"
                         "from which we want to recover contigs.")
parser.add_argument("-p", "--pos", required=True,
                    help="csv file of positions of contigs obtained with the"
                         "extract_pos_from_sam.py script [needs a .sam alignment input file,"
                         "obtained with e.g., BWA or GraphMap]")
parser.add_argument("-r", "--ref_folder", required=True,
                    help="path to directory that contains reference genome"
                         "in 17 separate .fsa files (see ref_names list in the code)."
                         "Download the files : http://www.yeastgenome.org/strain/S288C/overview")
parser.add_argument("-d", "--out_dir", default='./',
                    help="Directory for output and temporary files.")
parser.add_argument("-e", "--enzyme", default='BamHI',
                    help="Restriction enzyme to compute simulated optical map"
                         "for the yeast genome (chose in enz_dic dictionary in the code).")
parser.add_argument("--plot_length", type=bool, default=False,
                    help="whether to plot histogram of well/wrongly"
                         "ordered contigs in function of contig length"
                         "instead of number of restriction site, or not.")

args = parser.parse_args()

ref_fn = args.sequences
chr_pd_fn = args.pos
root_ref = args.ref_folder
root_out = args.out_dir
enz_name = args.enzyme
plot_length = args.plot_length

RE_SITE = enz_dic[enz_name]

ref_list = list(SeqIO.parse(ref_fn, "fasta"))
n_reads = len(ref_list)

minLen = 100

# Use synthetic optical map
ref_names = ["chr01.fsa", "chr02.fsa", "chr03.fsa", "chr04.fsa", "chr05.fsa",
             "chr06.fsa", "chr07.fsa", "chr08.fsa", "chr09.fsa",
             "chr10.fsa", "chr11.fsa", "chr12.fsa", "chr13.fsa",
             "chr14.fsa", "chr15.fsa", "chr16.fsa", "chrmt.fsa"]
numMaps = len(ref_names)
refMaps = {}
begin_ref_len = np.empty(numMaps)
k = 0
for ref_name in ref_names:
    seq_name = root_ref + ref_name
    fragment = SeqIO.parse(seq_name, "fasta").next().seq
    (refMaps[k], _, begin_ref_len[k], _) = make_fragment_maps(fragment, RE_SITE, minLen)
    k += 1

# Load reads' positions
whichChrTrue = np.zeros(n_reads)
positionTrue = np.zeros(n_reads)
num_chr = 17
alns_df = pd.DataFrame.from_csv(chr_pd_fn, sep='\t')
for ki in range(num_chr):
    chr_nb = int_to_roman(ki+1)
    sub_df = alns_df[alns_df['chr_nb'] == chr_nb]
    indexes = sub_df['number']
    whichChrTrue[indexes] = ki + 1
    sub_pos = sub_df['pos']
    positionTrue[indexes] = sub_pos


time1 = time.time()
sigmaInv = 1
C = 3000
bestPosDic = {}
bestPosA = np.empty(n_reads)
whichChr = np.empty(n_reads)
strandsA = np.empty(n_reads)
read_lengths = np.empty(n_reads)
out_fn = root_out + 'opt_map_%s.res' % (enz_name)
out_fh = open(out_fn, 'wb')
out_fh.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
('read nb', 'read id', 'bestMapChr', 'strand', 'bestPos', 'bestMapTrue', 'bestPosTrue'))
for readnb in range(n_reads):
    record = ref_list[readnb]
    read_lengths[readnb] = len(record.seq)
    scores = np.empty(2*numMaps)
    endPoss = np.empty(2*numMaps)
    read_id = record.id
    # read_id = ref.references[readnb]
    for mapNb in range(numMaps):
        refMap = refMaps[mapNb]
        fragment = record.seq
        (fragMap, fragMapRev, beginLength, endLength) = make_fragment_maps(fragment, RE_SITE, minLen)
        (endPos, score) = find_best_match_DP(fragMap, refMap, C, sigmaInv, endLength, 40)
        (endPosRev, scoreRev) = find_best_match_DP(fragMapRev, refMap, C, sigmaInv, beginLength, 40)

        scores[mapNb] = score
        endPoss[mapNb] = endPos
        scores[numMaps + mapNb] = scoreRev
        endPoss[numMaps + mapNb] = endPosRev

    bestMap = np.argmax(scores)
    bestPos = endPoss[bestMap]
    strand = +1
    if bestMap > numMaps -1:
        bestMap = bestMap - numMaps
        strand = -1
    # Add outer cut part of chromosome
    bestPos += begin_ref_len[bestMap]
    strandsA[readnb] = strand
    bestPosDic[readnb] = (bestMap, strand, bestPos)
    bestPosA[readnb] = bestPos
    whichChr[readnb] = bestMap +1
    # Compare to ``true'' posiiton and chromosome
    bestMapTrue = whichChrTrue[readnb]
    bestPosTrue = positionTrue[readnb]
    # print result to file
    out_fh.write('%d\t%s\t%d\t%d\t%d\t%d\t%d\n' % (readnb, read_id, bestMap, strand, bestPos, bestMapTrue, bestPosTrue))
out_fh.close()
print time.time() - time1

# Print results in a figure
lenMaps = np.empty(n_reads)
for readnb in range(n_reads):
    record = ref_list[readnb]
    fragment = record.seq
    (fragMap, fragMapRev, beginLength, endLength) = make_fragment_maps(fragment, RE_SITE, minLen)
    lenMaps[readnb] = len(fragMap)
    # Plot in function of contig length
    if plot_length:
        lenMaps[readnb] = 10*int(1e-4*len(record.seq))

idxs = np.argsort(-lenMaps)
chrSrt = whichChr[idxs]
chrTrueSrt = whichChrTrue[idxs]
posSrt = bestPosA[idxs]
posSrtTrue = positionTrue[idxs]


# Rewrite results sorted by number of occurrences
out_fn = root_out + 'opt_map_sorted_%s.res' % (enz_name)
out_fh = open(out_fn, 'wb')
out_fh.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('read nb.', 'read length', 'strand', 'nb. RS', 'Chr. (O.M.)', 'Chr. (BWA)', 'Pos. (O.M)', 'Pos (BWA)'))
for idx in range(n_reads):
    readnb = idxs[idx]
    # bestChrV = chrSrt[idx]
    bestChrV = whichChr[readnb]
    # trueChrV = chrTrueSrt[idx]
    trueChrV = whichChrTrue[readnb]
    # bestPosV = posSrt[idx]
    bestPosV = bestPosA[readnb]
    # truePosV = posSrtTrue[idx]
    truePosV = positionTrue[readnb]
    nb_RS_Occ = lenMaps[readnb]
    strand = strandsA[readnb]
    record = ref_list[readnb]
    read_len = len(record.seq)
    out_fh.write('%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n' % (readnb, read_len, strand, nb_RS_Occ, bestChrV, trueChrV, bestPosV, truePosV))
out_fh.close()

lenSort = np.empty(n_reads)
nbDifs = np.empty(n_reads)
for i in range(n_reads):
    nbDifs[i] = np.count_nonzero(np.array(chrSrt[:i+1]) - np.array(chrTrueSrt[:i+1]))
    lenSort[i] = lenMaps[idxs[i]]
    # For contig lengt instead of number of RS
    # lenSort[i] = int(0.001*read_lengths[idxs[i]])

nbDif = np.count_nonzero(chrSrt - chrTrueSrt)

lenL = np.empty(len(set(lenSort)))
trueL = np.empty(len(set(lenSort)))
falseL = np.empty(len(set(lenSort)))
cpt = 0
for num in set(lenSort):
    ind = np.argwhere(lenSort == num)
    ind = [int(el) for el in ind]
    lenL[cpt] = num
    falseL[cpt] = np.count_nonzero(np.array(chrSrt[min(ind):max(ind)+1]) - np.array(chrTrueSrt[min(ind):max(ind)+1]))
    trueL[cpt] = len(ind) - falseL[cpt]
    cpt += 1

# width = 0.35
width = 2./100*max(lenL[2:])

plt.ylabel('number of contigs', fontsize=18)
plt.xlabel('number of RS', fontsize=18)
namefig = root_out + 'OpticalMapSimulatedS288C_' + enz_name + '.eps'

if plot_length:
    plt.xlabel('contig length (kb)', fontsize=18)
    namefig = root_out + 'OpticalMapSimulatedS288C_' + enz_name + '_len.eps'

p1 = plt.bar(lenL[2:], trueL[2:], width, color=(0.2588,0.4433,1.0))#, edgecolor='none')
p2 = plt.bar(lenL[2:], falseL[2:], width, color=(1.0,0.5,0.62), bottom=trueL[2:])#, edgecolor='none')

plt.legend((p1[0], p2[0]), ('correctly mapped', 'mis-mapped'))
plt.savefig(namefig, bbox_inches='tight')

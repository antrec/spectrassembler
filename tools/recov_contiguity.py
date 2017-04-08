# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 11:58:58 2015

@author: antlaplante
"""
from Bio import SeqIO
import subprocess
import argparse



def write_ends(reads_fn, ends_fn, end_len):

    if reads_fn[-1] == 'q':
        fmt = 'fastq'
    else:
        fmt = 'fasta'

    len_dic = {}
    fh = open(reads_fn, 'rU')
    ends_fh = open(ends_fn, 'wb')

    for record in SeqIO.parse(fh, fmt):
        read_id = record.id
        seq = record.seq
        len_dic[read_id] = len(seq)
        trim_len = min(end_len, len(seq))
        begin_seq = seq[:trim_len]
        end_seq = seq[len(seq)-trim_len:]

        ends_fh.write(">%s\n%s\n" % (read_id + "_b", begin_seq))
        ends_fh.write(">%s\n%s\n" % (read_id + "_e", end_seq))

    fh.close()
    ends_fh.close()

    return len_dic

def run_minimap(minimap_path, seqs_fn, out_fn):

    out_fh = open(out_fn, 'wb')

    cmd = [minimap_path, '-Sw5', '-L', '500', seqs_fn, seqs_fn]

    subprocess.call(cmd, stdout=out_fh)

    out_fh.close()

    return

def retrieve_ovl_fn(ovl_fn, ends_ovl_fn, len_dic, end_len, max_hang):
    temp_fh = open(ends_ovl_fn, "r")
    ovl_fh = open(ovl_fn, "wb")
    for line in temp_fh:
        fields = line.split()
        id1 = fields[0][:-2]
        id2 = fields[5][:-2]
        strand = fields[4]
        p1 = fields[0][-1]
        p2 = fields[5][-1]
        if id1 == id2:
            continue
        astart = int(fields[2])
        aend = int(fields[3])
        bstart = int(fields[7])
        bend = int(fields[8])
        alen = len_dic[id1]
        blen = len_dic[id2]

        # Added exclusion of overlaps that are not really in the ends
        a_begin = (astart < max_hang)
        a_end = (aend > end_len - max_hang)
        b_begin = (bstart < max_hang)
        b_end = (bend > end_len - max_hang)

        same_ends = (a_begin and b_begin) or (a_end and b_end)
        opp_ends = (a_begin and b_end) or (a_end and b_begin)

        is_overlap = ((strand == '+') and (p1 != p2) and opp_ends) or \
                     ((strand == '-') and (p1 == p2) and same_ends)

        if not is_overlap:
            continue

        if p1 == 'e':
            astart = min(alen, alen + astart - end_len)
            aend = min(alen, alen + aend - end_len)
        if p2 == 'e':
            bstart = min(blen, blen + bstart - end_len)
            bend = min(blen, blen + bend - end_len)

        ovl_fh.write("%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n" %
                    (id1, alen, astart, aend, fields[4], id2, blen, bstart, bend, fields[9], fields[10], fields[11], fields[12]))

    ovl_fh.close()
    temp_fh.close()
    return

def run_miniasm(miniasm_path, reads_fn, ovl_fn, assembly_gfa):

    assembly_fh = open(assembly_gfa, 'wb')

    cmd = [miniasm_path, '-e', '0', '-1', '-2', '-c', '0', '-r', '1,0.', '-f', reads_fn, ovl_fn]

    subprocess.call(cmd, stdout=assembly_fh)

    assembly_fh.close()

    return

def get_new_contigs(gfa_fn, reads_fn, new_contigs_fn):

    with open(new_contigs_fn, 'wb') as new_fh:
        cmd = ['awk', '/^S/{print \">\"$2\"\\n\"$3}', gfa_fn]
        awk = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        subprocess.call(['fold'], stdin=awk.stdout, stdout=new_fh)
        awk.wait()

    included_seqs = []
    gfa_fh = open(gfa_fn, 'rb')
    for line in gfa_fh.readlines():
        if line[0] == 'a':
            included_seqs.append(line.split()[3])

    gfa_fh.close()

    if reads_fn[-1] == 'q':
        fmt = 'fastq'
    else:
        fmt = 'fasta'

    with open(new_contigs_fn, 'a') as new_fh:
        with open(reads_fn, 'rU') as reads_fh:
            for record in SeqIO.parse(reads_fh, fmt):
                if not(record.id in included_seqs):
                    new_fh.write(record.format(fmt))

    return


parser = argparse.ArgumentParser(description="Uses miniasm to stick together"
                                             "contigs from spectrassembler."
                                             "Needs minimap and miniasm.")
parser.add_argument("sequences",
                    help="sequences output from spectrassembler in FASTA format"
                         "from which we want to recover contigs.")
parser.add_argument("--minimap", default="minimap",
                    help="path/to/minimap")
parser.add_argument("--miniasm", default="miniasm",
                    help="path/to/miniasm")
parser.add_argument("-e", "--end_len", type=int, default=40000,
                    help="length of ends of large sequences"
                         "to compare for overlaps.")
parser.add_argument("-h", "--max_hang", type=int, default=300,
                    help="max distance to end of sequence"
                         "for an overlap to be kept.")

args = parser.parse_args()

minimap_path = args.minimap_path
miniasm_path = args.miniasm_path
end_len = args.end_len
max_hang = args.max_hang
reads_fn = args.sequences
name = reads_fn.split('.')[:-1]

ends_fn = '%s.ends.%s' % (name, reads_fn.split('.')[-1])
len_dic = write_ends(reads_fn, ends_fn, end_len)

ends_ovl_fn = '%s.ends.paf' % name
run_minimap(minimap_path, ends_fn, ends_ovl_fn)

ovl_fn = '%s.fromends.paf' % name
retrieve_ovl_fn(ovl_fn, ends_ovl_fn, len_dic, end_len, max_hang)

assembly_fn = '%s.fromends.gfa' % name
run_miniasm(miniasm_path, reads_fn, ovl_fn, assembly_fn)

glued_contigs_fn = '%s.contigs.fasta' % name
get_new_contigs(assembly_fn, reads_fn, glued_contigs_fn)

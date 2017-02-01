#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""
Check that there are no repeated read id in the data.
If there are, create a new data file [called reads.pp.fasta if the
old one was reads.fasta for example] where the reads ids are
modified in order to have no id duplicates.
"""
from Bio import SeqIO
import sys


def change_name(record, new_name):
    """
    Changes the id, name and description of a SeqIO record.

    Parameters
    ----------
    record : SeqIO record
    new_name : str

    """
    record.id = new_name
    record.name = new_name
    descript_list = record.description.split(' ')
    descript_list[0] = new_name
    record.description = ' '.join(descript_list)

# Load reads
fn = sys.argv[1]
fh = open(fn, "rU")
fmt = fn.split('.')[-1]
print u'open reads file {} to check that read ids are unique'.format(fn)
record_list = list(SeqIO.parse(fh, fmt))
fh.close()

# Construct {read name : read number} dictionary
read_nb_dic = {}
cpt = 0
nb_doublon = {}
for (idx, record) in enumerate(record_list):
    if read_nb_dic.has_key(record.id):
        print u'id {} same for reads {} and {}'.format(record.id, read_nb_dic[record.id], cpt)
        if nb_doublon.has_key(record.id):
            nb_doublon[record.id] += 1
        else:
            nb_doublon[record.id] = 1
        new_id = record.id + '_%d' % nb_doublon[record.id]
        change_name(record_list[idx], new_id)

    read_nb_dic[record_list[idx].id] = cpt
    cpt += 1

if len(nb_doublon.keys()) > 0:
    out_fn = '.'.join(fn.split('.')[:-1]) + '.pp.' + fmt
    print u'write renamed reads with unique ids to file {}'.format(out_fn)
    with open(out_fn, "w") as handle:
        SeqIO.write(record_list, handle, fmt)
else:
    print u'Reads ids are already unique in file {}. No changes to do.'.format(fn)

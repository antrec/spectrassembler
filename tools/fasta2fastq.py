#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""
Create a fastq file from fasta file with fake quality values all equal.
"""
import sys
from Bio import SeqIO


# Get inputs
fa_path = sys.argv[1]
fq_path = sys.argv[2]

# Make fastq
with open(fa_path, "rb") as fasta, open(fq_path, "wb") as fastq:
    for record in SeqIO.parse(fasta, "fasta"):
        record.letter_annotations["phred_quality"] = [40] * len(record)
        SeqIO.write(record, fastq, "fastq")

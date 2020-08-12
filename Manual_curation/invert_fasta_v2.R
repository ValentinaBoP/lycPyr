#!/usr/bin/python

import sys
from Bio.Seq import Seq
from pyfaidx import Fasta

prefix = str(sys.argv[1])
input_fasta = prefix + ".fasta"
input_coord = prefix + ".coordinates"
with open (input_coord) as coords:
        with Fasta(input_fasta, mutable = True) as fasta:
                for record in fasta:
                        for line in coords:
                                begin, end = line.rstrip().split()[1:3]
                                begin = int(begin)
                                end = int(end)
                                seq = record[begin:end]
                                #inv = str(seq[::-1])
                                #inv = ''.join(inv)
                                inv = str(Seq(str(seq)).reverse_complement())
                                record[begin:end] = inv

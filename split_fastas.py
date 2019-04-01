#!/usr/bin/python2.7
#RANDOMLY split reads based on a "split" factor as stated in the call as the third variable, see line 26. also see subsample_random.py
#call as (e.g. to split in two): python split_fastas.py all_reads.fasta outfile_A.fasta outfile_B.fasta 2

import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import random

#handle input arguments
import sys

infile = sys.argv[1]
outfile_A = sys.argv[2]
outfile_B = sys.argv[2]
splits = sys.argv[3]


#read in fasta and prep number of necessary outfiles
record = list(SeqIO.parse(infile, "fasta"))
outA = list()
outB = list()

#prep the random indexing and size of subsample
x = (range(0, len(record)))
random.shuffle(x)
y = len(record)*(1/splits)
z = int(round(y))
splitA = x[0:z]
splitB = x[z:len(record)]

outA = record[splitA]
outB = record[splitB]

#safe as fastas again
SeqIO.write(outA, outfile_A, "fasta")
SeqIO.write(outB, outfile_B, "fasta")
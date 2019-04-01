#!/usr/bin/python2.7
#RANDOMLY subsample sequencing files to e.g. half their number of reads, set factor in line 25
#call as (e.g. to subset half "0.5" the reads): python2.7 subsample_random.py all_reads.fasta subset_reads.fasta #0.5

import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import random

#handle input arguments
import sys

infile = sys.argv[1]
outfile = sys.argv[2]
#factor = sys.argv[2]


# in a loop
record = list(SeqIO.parse(infile, "fasta"))
out = list()

#prep the random indexing and size of subsample
x = (range(0, len(record)))
random.shuffle(x)
y = len(record)*0.5 #randomisation factor
z = int(round(y))
subset = x[0:z]

#extract subset in a loop of randomly picked index numbers to defined size/ factor
for n in subset: 
 out.append(SeqRecord(record[n].seq, id = record[n].id))

#safe as a fasta again
SeqIO.write(out, outfile, "fasta")

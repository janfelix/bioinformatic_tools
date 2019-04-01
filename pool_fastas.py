#!/usr/bin/python2.7
#pool sequences from multiple fasta files (using * wild card) into one fasta file, easily incorporate a processing step in line 18...
#call as pool_fastas.py outfile.fasta

#handle input arguments
import sys
outfile_1 = sys.argv[1]

import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
record1 = list()

import glob
for filename in glob.glob('*_pool.fasta'):
	reads = list(SeqIO.parse(filename, "fasta"))
	for n in xrange(0, len(reads)):
		record1.append(SeqRecord(reads[n].seq, id = reads[n].id))

SeqIO.write(record1, outfile_1, "fasta")
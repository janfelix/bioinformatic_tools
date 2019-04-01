#!/usr/bin/python2.7
#subset sequences such as OTUs based on a list of labels given i
#call as: python2.7 subset_fastas.py inputfile.fasta wantedreads.csv goodreads.fasta badreads.fasta



#handle input arguments, fancy option for later
import sys
infile_1 = sys.argv[1]
infile_2 = sys.argv[2]
outfile_1 = sys.argv[3]
outfile_2 = sys.argv[4]


#open necessary modules
import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re
import csv


#open sequences fasta as a list
record = list(SeqIO.parse(infile_1, "fasta"))

#open hits txt with as a list of lists
with open(infile2, 'r') as myfile:
	wr = csv.reader(myfile, delimiter='\t')
	otus = list(wr)
#extract first item from each sub-list, the labels and save as separate list
labels = [item[0] for item in otus]

#open empty lists to save to
out1 = list()
out2 = list()

#loop through each sequence, check against the list of hit labels
for n in xrange(0, len(record)):
	if record[n].id in labels: out1.append(SeqRecord(record[n].seq, id = record[n].id))
	else: out2.append(SeqRecord(record[n].seq, id = record[n].id))

#save output as fasta, save as goodreads when matching and badreads when not
#SeqIO.write(out1, "goodreads.fasta", "fasta")
SeqIO.write(out1, outfile_1, "fasta")

#SeqIO.write(out2, "badreads.fasta", "fasta")
SeqIO.write(out2, outfile_2, "fasta")

#!/usr/bin/python3
#extract wanted sequences from a fasta file based on a list of defined sequences in a csv file eg from a blast result
#call as: python extract_hits.py infile1.csv infile2.fasta outfile1.fasta outfile2.fasta

import sys
#handle input arguments
infile1 = sys.argv[1]
infile2 = sys.argv[2]
outfile1 = sys.argv[3]
outfile2 = sys.argv[4]

import pandas
import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#get list of nodes with hits in diamond search
hits = pandas.read_csv(infile1)
hl = list(hits['qseqid'])

#open empty lists to save to
out1 = list()
out2 = list()

#open transcript fasta file from spades assembler
record = list(SeqIO.parse(infile2, "fasta"))

#loop through each sequence, check against the list of hit labels; python automatically only extracts unique labels matching the number of set(hl)
for n in range(0, len(record)):
	if record[n].id in hl: out1.append(SeqRecord(record[n].seq, id = record[n].id))
	else: out2.append(SeqRecord(record[n].seq, id = record[n].id))
	
#save output as fasta, save wanted reads when matching in outfile1 and unwanted reads in outfile2
SeqIO.write(out1, outfile1, "fasta")
SeqIO.write(out2, outfile2, "fasta")






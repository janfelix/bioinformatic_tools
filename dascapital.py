#!/usr/bin/python2.7
#Capitalise sequences in fasta files (if necessary for down stream analysis), and capitalize on it
#call as: python dascapital.py infile.fasta outfile.fasta


#handle input arguments
import sys

infile_1 = sys.argv[1]
outfile_1 = sys.argv[2]

#read infile
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

record1 = list(SeqIO.parse(infile_1, "fasta"))
record3 =list()

#getting err done
for n in range(0, len(record1)):record3.append(SeqRecord(Seq(str(record1[n].seq).upper()), id=record1[n].id))



#safe as a fasta again
SeqIO.write(record3, outfile_1, "fasta")
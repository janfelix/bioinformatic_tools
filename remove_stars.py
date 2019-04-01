#!/usr/bin/python2.7
#just remove stars from aa sequences and add a sample_id
#call as: python2.7 remove_stars.py translated.fasta sample_id


#handle input arguments
import sys

infile_1 = sys.argv[1]
sample_id = sys.argv[2]

#read infile
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

record1 = list(SeqIO.parse(infile_1, "fasta"))
record3 =list()

#getting err done
for n in xrange(0, len(record1)):record3.append(SeqRecord(Seq(str(record1[n].seq).replace("*","")), id=sample_id+'_'+record1[n].id))


#safe as a fasta again
SeqIO.write(record3, "aa_no_*.faa", "fasta")
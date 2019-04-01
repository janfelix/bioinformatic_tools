#!/usr/bin/python2.7
#read-in aa fasta files, select for min length and remove eg stars (line 24), add a sample id,  save as fasta output for further processing
#call as: python2.7 clean_fasta.py translated.fasta sample_id


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
for n in xrange(0, len(record1)):
	if len(record1[n].seq) >= 140: record3.append(SeqRecord(Seq(str(record1[n].seq).replace("*","")), id=sample_id+'_'+record1[n].id))

	

#safe as a fasta again
SeqIO.write(record3, "aa_clean.fasta", "fasta")
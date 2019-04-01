#!/usr/bin/python2.7
#call as: python2.7 cluster_size_overview.py clusterfile.fasta cluster_labels.csv

#handle input arguments
import sys

infile_1 = sys.argv[1]
outfile_1 = sys.argv[2]

import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re

record1 = list(SeqIO.parse(infile_1, "fasta"))

import csv
#prep csv outputfile
# first we write the header :
header = [['cluster', 'label', 'size']]

#create csv
with open(outfile_1, 'w') as outcsv:
 a = csv.writer(outcsv)
 a.writerows(header)

#loop through the file in the handle
 for p in range(0, len(record1)):
 
        cluster = record1[p].id
        label, slice, size, rice = re.split('=|;', record1[p].id)
        
        # build the list for the cluster
        list_cluster = [[cluster, label, size]]
        #write to csv file
        a.writerows(list_cluster)






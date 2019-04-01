#!/usr/bin/python2.7
#extract statistics on cluster sizes from sequencing clusters
#call as: python2.7 cluster_size_overview.py clusterfile.fasta clusterstats.csv

#handle input arguments
import sys

infile_1 = sys.argv[1]
outfile_1 = sys.argv[2]

import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re

# in a loop
record = list(SeqIO.parse(infile_1, "fasta"))
out = list()

for n in xrange(0, len(record)):
	crap, slice, size, rice = re.split('=|;', record[n].id)
	out.append(int(size))

import numpy

size = numpy.array([out])

_0 = [0,len(out),sum(out)]
_1 = [1,len(size[size>1]),sum(size[size>1])]
_10 = [10, len(size[size>10]), sum(size[size>10])]
_100 = [100, len(size[size>100]), sum(size[size>100])]
_500 = [500,len(size[size>500]),sum(size[size>500])]
_1000 = [1000, len(size[size>1000]), sum(size[size>1000])]
_10000 = [10000, len(size[size>10000]), sum(size[size>10000])]
_100000 = [100000, len(size[size>100000]), sum(size[size>100000])]

#result = [_0, _1, _10, _100, _1000, _10000]

result = numpy.asarray([_0, _1, _10, _100, _500, _1000, _10000, _100000])
numpy.savetxt(outfile_1, result, delimiter=",")


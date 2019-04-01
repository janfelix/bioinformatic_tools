#!/usr/bin/python3
#convert blast and diamond blast output txt files to csv files
#call as: python tab2csv.py infile.txt outfile.csv

import sys
import pandas as pd

#handle input arguments
infile = sys.argv[1]
outfile = sys.argv[2]
#define header
cols = ["qseqid", "qlen", "qstart", "qend", "sseqid", "sstart", "send", "length", "pident", "mismatch", "gapopen", "evalue", "bitscore", "score", "staxids", "stitle"]
#load txt
#df = pd.read_csv(infile, sep='\s+', header=None, names=cols)
df = pd.read_csv(infile, sep='\t', header=None, names=cols)
#save as csv
df.to_csv(outfile)

#!/usr/bin/python3
#extract blast hits with the highest pident (line 27) for a given taxon for each query sequence
#call as: python taxon_maxid_csv.py infile.csv taxon outfile.csv

import sys
import pandas as pd
import re

#handle input arguments
infile = sys.argv[1]
taxon = sys.argv[2]
outfile = sys.argv[3]

#load file and extract taxon hits
df = pd.read_csv(infile)
dfk = df[df['stitle'].str.contains(taxon, case=False, regex=False)]
#open empty df
dfo = pd.DataFrame()

#scan hits for top hit per taxon
qseqs = dfk['qseqid'].unique().tolist()

#extract and split data, 'qmax.iloc' used instead of eg 'str(qmax['qseqid'])' to get the full length strings from df cells even if very long
for q in qseqs:
 
        dfq = dfk[dfk['qseqid']==q]
        qmax = dfq[dfq['pident']==dfq['pident'].max()]
        NODE, nodenum, length, lennum, cov, covnum, g, i = re.split('_', qmax.iloc[0, 1])
        qmax = qmax.assign(qnode=nodenum)
        qmax = qmax.assign(qcove=covnum)
        taxn = re.findall(r'\[(.*?)\]', qmax.iloc[0, 16])
        if len(taxn)==0:
        	taxn ="none"
        else: 
        	taxn=taxn[0]
        qmax = qmax.assign(stax=taxn)
        dfo = dfo.append(qmax)

#add sample column 
sample, scheiss = re.split('_', infile)
dfo['sample']= sample
   
#save as csv
dfo.to_csv(outfile)

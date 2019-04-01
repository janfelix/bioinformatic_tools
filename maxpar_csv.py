#!/usr/bin/python3
#extract blast or diamond blast results for sequences for the top hit based on eg the bitscore as defined by "parameter". Optionally subset data based on a search phrase, see line 16
#call as: python maxpar_csv.py infile.csv parameter outfile.csv

import sys
import pandas as pd
import re

#handle input arguments
infile = sys.argv[1]
par = sys.argv[2]
outfile = sys.argv[3]

#load file, in this version keep all taxon hits
df = pd.read_csv(infile)
#dfk = df[df['stitle'].str.contains('phrase|alternative_phrase', case=False, regex=True)]
dfk = df
#open empty df
dfo = pd.DataFrame()

#scan hits for top hit per taxon
qseqs = dfk['qseqid'].unique().tolist()

#extract and split data, 'qmax.iloc' used instead of eg 'str(qmax['qseqid'])' to get the full length strings from df cells even if very long
for q in qseqs:
 
        dfq = dfk[dfk['qseqid']==q]
        qmax = dfq[dfq[par]==dfq[par].max()]
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

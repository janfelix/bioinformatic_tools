#!/usr/bin/python3
#take highest match by bitscore from blast analyses and extract the taxonomic lineage from the ncbi taxdump database 
#call as: python max2lineage_multi.py infile.txt outfile.csv numcores
import pandas as pd
import re
import sys
import multiprocessing
from multiprocessing import Pool

#handle input arguments, define lineage db
infile = sys.argv[1]
infile1 = '/data3/ncbi_taxdump_20220503/rankedlineage.dmp'
outfile = sys.argv[2]
pools = sys.argv[3]

#load txt, define header and format staxids column
cols = ["qseqid", "qlen", "qstart", "qend", "sseqid", "sstart", "send", "length", "pident", "mismatch", "gapopen", "evalue", "bitscore", "score", "staxids", "sscinames", "skingdoms", "sphylums", "stitle"]
df = pd.read_csv(infile, sep='\t', header=None, names=cols)
if df.staxids.dtype=='float64':
    df['staxids'] = df['staxids'].fillna(float(1))
    df['staxids'] = df['staxids'].astype(int)
    df['staxids'] = df['staxids'].astype(str)
else:
    df['staxids'] = df['staxids'].fillna(str(1))
print("infile loaded")

#read lineage db into list of lists >> using pandas, use only columns with actual data and not "|" when building the list
df1= pd.read_csv(infile1, sep='\t|\t', header=None)
taxlin=df1.iloc[:,[0,2,4,6,8,10,12,14,16,18]].values.tolist()
print("lineage db loaded")

#extract querry ids
qseqs = df['qseqid'].unique().tolist()
print("number of qseqs: "+str(len(qseqs)))

###multiprocessing
#build function
def lineage(q):
    dfq = df[df['qseqid']==q]
    qmax = dfq[dfq['bitscore']==dfq['bitscore'].max()]
    j = re.split(';', str(qmax['staxids'].values[0]))[0]
    lineage=[[id,tx,sp,ge,fa,od,cl,ph,ki,sk] for (id,tx,sp,ge,fa,od,cl,ph,ki,sk) in taxlin if int(id) == int(j)]
    if len(lineage) ==0:
        lineage = ['id','tx','sp','ge','fa','od','cl','ph','ki','sk']
    else:
        lineage = lineage[0] 
    return qmax.iloc[0].tolist()+lineage

#pool processes in parallel
po=Pool(int(pools))
data=po.map(lineage,qseqs)
po.close()

###combine to df, save as csv
df2 = pd.DataFrame(data, columns = ["qseqid", "qlen", "qstart", "qend", "sseqid", "sstart", "send", "length", "pident", "mismatch", "gapopen", "evalue", "bitscore", "score", "staxids", "sscinames", "skingdoms", "sphylums", "stitle", "stax_id", "tax_name", "species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"])
#df2.columns = ["qseqid", "qlen", "qstart", "qend", "sseqid", "sstart", "send", "length", "pident", "mismatch", "gapopen", "evalue", "bitscore", "score", "staxids", "sscinames", "skingdoms", "sphylums", "stitle", "stax_id", "tax_name", "species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"]
df2.to_csv(outfile)
print("done!")
#!/usr/bin/python3
#Takes the max hit per query based on the bitscore and retrieves the taxonomic lineage from NCBI's taxdump db
#call as: python max2lineagex.py infile.txt lineage_db outfile.csv
import pandas as pd
import re
import sys

#handle input arguments, define lineage db matching diamond blast db
infile = sys.argv[1]
infile1 = sys.argv[2]
outfile = sys.argv[3]

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

#open empty list
dfx = list()

#scan hits for top hit per taxon and extract lineage based on staxid
qseqs = df['qseqid'].unique().tolist()
print("number of qseqs: "+str(len(qseqs)))

for q in qseqs:
    dfq = df[df['qseqid']==q]
    qmax = dfq[dfq['bitscore']==dfq['bitscore'].max()]
    j = re.split(';', str(qmax['staxids'].values[0]))[0]
    lineage=[[id,tx,sp,ge,fa,od,cl,ph,ki,sk] for (id,tx,sp,ge,fa,od,cl,ph,ki,sk) in taxlin if int(id) == int(j)]
    if len(lineage) ==0:
        lineage = ['id','tx','sp','ge','fa','od','cl','ph','ki','sk']
    else:
        lineage = lineage[0] 
    dfx.append(qmax.iloc[0].tolist()+lineage)

dfo = pd.DataFrame(dfx, columns = ["qseqid", "qlen", "qstart", "qend", "sseqid", "sstart", "send", "length", "pident", "mismatch", "gapopen", "evalue", "bitscore", "score", "staxids", "sscinames", "skingdoms", "sphylums", "stitle", "stax_id", "tax_name", "species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"])
dfo.to_csv(outfile)
print("done!")
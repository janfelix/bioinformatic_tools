###call as: python taxid2lineage.py rankedlineage.dmp infile.csv outfile.csv

import pandas as pd
import re
import sys


infile1 = sys.argv[1]
infile2 = sys.argv[2]
outfile1 = sys.argv[3]



###read lineage db into list of lists>>delimiter problem!!!
#using pandas, use only columns with actual data and not the stupid "|" when building the list
df1= pd.read_csv(infile1, sep='\t|\t', header=None)
taxlin=df1.iloc[:,[0,2,4,6,8,10,12,14,16,18]].values.tolist()

###extract transcriptome annotation from csv
#load txt with diamond hits of contigs, fill nan values with 1, the lineage root
#df2 = pd.read_csv(infile2, sep=',', header=None, names=cols)
df2 = pd.read_csv(infile2, sep=',', header=0)
df2['staxids'] = df2['staxids'].fillna(str(1))


###extract staxid for one sample or multiple sample depending on the matrix structure
staxid = df2['staxids'].tolist()

###open empty data frame for lineages 
dfo = pd.DataFrame()


###extract lineage for staxids
#append lineage to taxon data
for i in staxid:

        j = re.split(';', i)[0]
        lineage=[[id,tx,sp,ge,fa,od,cl,ph,ki,sk] for (id,tx,sp,ge,fa,od,cl,ph,ki,sk) in taxlin if int(id) == int(j)]
        if len(lineage) ==0:
        	    lineage = ['id','tx','sp','ge','fa','od','cl','ph','ki','sk']
        else:
        	    lineage = lineage[0] 
        dfo = dfo.append(pd.DataFrame(data=[lineage]))

#define additional column header
dfo.columns = ["stax_id", "tax_name", "species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"]

###combine df4 and dfo, save as csv
df2.reset_index(drop=True, inplace=True)
dfo.reset_index(drop=True, inplace=True)
df5=pd.concat([df2, dfo], axis=1)
df5.to_csv(outfile1)
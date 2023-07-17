#!/usr/bin/env python3
#Translates transcripts into six frame amino acid sequences without the start and stop codons for partial metatranscriptomes
#call as:python3 66.py -n nucl.fasta -t 1 -a aa.fasta 
#translation_table "-t" corresponds to genetic codes NCBI, defaults to "1" (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)

import os
import optparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

inargs = optparse.OptionParser()
inargs.add_option("-n",  action = "store", type = "string", dest = "nucl")
inargs.add_option("-t",  action = "store", type = int, default=1, dest = "table")
inargs.add_option("-a",  action = "store", type = "string", dest = "aa")
(options, args) = inargs.parse_args()
record = list(SeqIO.parse(options.nucl, "fasta"))
ttable=int(options.table)

def six(record, ttable):
    fasta.append(SeqRecord(record.seq[:].ungap('-').translate(table=ttable),id=record.id+'_'+str('f0')))
    fasta.append(SeqRecord(record.seq[1:].ungap('-').translate(table=ttable),id=record.id+'_'+str('f1')))
    fasta.append(SeqRecord(record.seq[2:].ungap('-').translate(table=ttable),id=record.id+'_'+str('f2')))
    fasta.append(SeqRecord(record.seq[:].ungap('-').reverse_complement().translate(table=ttable),id=record.id+'_'+str('r0')))
    fasta.append(SeqRecord(record.seq[:-1].ungap('-').reverse_complement().translate(table=ttable),id=record.id+'_'+str('r1')))
    fasta.append(SeqRecord(record.seq[:-2].ungap('-').reverse_complement().translate(table=ttable),id=record.id+'_'+str('r2')))

fasta=list()
for n in range(0, len(record)):
    six(record[n],ttable)

SeqIO.write(fasta, options.aa, "fasta")
#merge genome and protein CDS fasta files into a genbank file, some editing in the outputfile might be required
#call as: python fasta2genbank.py dnaseqfile.fasta cdsfile.fasta genbankfile.gb
#cds locations must be in the name following a ..._range_start_stop pattern(see line 33)
#feature description must be in the name within <...> (see line34)

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna, generic_protein
from Bio.SeqFeature import SeqFeature, FeatureLocation
import re

#handle input arguments
seq_file = sys.argv[1]
cds_file= sys.argv[2]
genbank_file = sys.argv[3]

#read in DNA origin sequence
sequence = list(SeqIO.parse(seq_file, "fasta"))
for seq in sequence:seq.seq.alphabet = generic_dna

#open a record and add DNA sequence
sequence_object = Seq(str(sequence[0].seq), IUPAC.ambiguous_dna)
record = SeqRecord(sequence_object, id=sequence[0].id, name=sequence[0].name, description=sequence[0].description)

#read in protein CDS features
cds = list(SeqIO.parse(cds_file, "fasta"))
for seq in cds:seq.seq.alphabet = generic_protein 

#add features to the record              
for feat in cds: 
	n, nodenum, r, start, stop = re.split('_', feat.name)
	prod = re.findall(r'\<(.*?)\>', feat.description)
	feature = SeqFeature(FeatureLocation(start=int(start), end=int(stop)), type='CDS') #edit feature type manually in genbank file if necessary
	feature.qualifiers = {'product':prod, 'translation':str(feat.seq)}
	record.features.append(feature)

SeqIO.write(record, genbank_file, "genbank")


#call as: python orf_finder.py infile.fasta translation_table min_protein_length outfile.csv fasta_outfile.fasta
#translation_table corresponds to genetic codes at NCBI (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
#example: python orf_finder.py metagenome.fasta 1 100 metagenome_orf.csv metagenome_orf.fasta

#required packages
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import sys
import pandas
import re

#input variables
infile = list(SeqIO.parse(sys.argv[1], "fasta"))
table = int(sys.argv[2])
min_pro_len = int(sys.argv[3])
outfile = sys.argv[4]
fastafile = sys.argv[5]
#fasta list to store orfs
fasta = list()

#define function
def find_orfs_with_trans(record, trans_table, min_protein_length):
    answer = []
    seq = record.seq
    query = record.description
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end-aa_start >= min_protein_length:
                    if strand == 1:
                        start = frame+aa_start*3
                        end = min(seq_len,frame+aa_end*3+3)
                    else:
                        start = seq_len-frame-aa_end*3-3
                        end = seq_len-frame-aa_start*3
                    answer.append((query, start, end, strand, len(trans[aa_start:aa_end]), trans[aa_start:aa_end]))
                    fasta.append(SeqRecord(Seq(trans[aa_start:aa_end], IUPAC.protein), id = re.split(' ', query)[0]+'_'+'range'+'_'+str(start)+'_'+str(end)))
                aa_start = aa_end+1
    answer.sort()
    return answer

#run the function(s), in a loop feeding into a pandas data frame
df = pandas.DataFrame()
for p in infile:
	orf_list = find_orfs_with_trans(p, table, min_pro_len)
	df = df.append(orf_list)

#save data frame as csv
results = pandas.DataFrame(orf_list, columns=['query', 'start', 'end', 'strand', 'length', 'trans'])
df.to_csv(outfile, header=['query', 'start', 'end', 'strand', 'length', 'trans'])

#save fasta file with orfs
SeqIO.write(fasta, fastafile, "fasta")

print("YOU ORF!!!")




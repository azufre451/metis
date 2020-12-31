#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("fastaFile", help="FASTA file containing the sequences")
parser.add_argument("metadataFile", help="Metadata file containing the info")
parser.add_argument("--type", default="fasta", help="input file type (fasta|fastq)") 
parser.add_argument("--idField",default=0, type=int)
parser.add_argument("--separator",default='|', help="separator on SeqID for metadata matching")
parser.add_argument("--minseqs",default=0, type=int, help="Minimum Number of Sequences in the file")
parser.add_argument("--nometa",default='unknown', help="Label for samples that do not appear in metadata")
args=parser.parse_args()

ww={}
for line in open(args.metadataFile):
	ww[line.strip().split('\t')[0]] = line.strip().split('\t')[1::]


w={'A':'0001','T':'0010','C':'0100','G':'1000','-':'0000'}

rows=[]

for seq_record in SeqIO.parse(args.fastaFile, args.type):
	strr=[] 
	for chara in str(seq_record.seq):
		if chara.upper() in w.keys(): strr += list(w[chara.upper()])
		else: strr += list('0000')
	if args.separator: tgt=seq_record.id.split(args.separator)[0] 
	else: tgt = seq_record.id

	if tgt in ww:
		l= [ww[tgt][args.idField]]+ strr
	else:
		if args.nometa == 'SAMPLENAME': l= [seq_record.id]+ strr
		else: l =[args.nometa] + strr

	
	rows.append(l)

if len(rows) > args.minseqs: 

	prepender = ['#meta_1']+['g_'+str(x) for x in range(0,len(rows[0]))] 

	rows.insert(0,prepender)


	for nrow in [list(x) for x in zip(*rows)]:
		print '\t'.join(nrow)


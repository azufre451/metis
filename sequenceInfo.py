from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import argparse
import bz2,sys

parser = argparse.ArgumentParser(description="return stats for fastx sequences (Percent. of each nucleotide, length)")
parser.add_argument("fastaFile", help="FASTA file containing the sequences")
parser.add_argument("--type", default="fasta", help="input file type (fasta|fastq)")
parser.add_argument("--verbose", action="store_true")
parser.add_argument("--tabular", action="store_true")
parser.add_argument("--trimlen")
parser.add_argument("--snpsFile")
parser.add_argument("--bz2", action="store_true")
args=parser.parse_args()
 
seqList = []
total = 0
avn = 0
Nc = 0
li={'A':0,'T':0,'C':0,'G':0,'N':0,'-':0,'*':0}

if args.bz2: st=bz2.BZ2File(args.fastaFile,'r')
elif args.fastaFile == '-': st=sys.stdin
else: st=open(args.fastaFile,'r')

for seq_record in SeqIO.parse(st, args.type):
	total+=len(seq_record.seq)
	for i in seq_record.seq:
		if i.upper() in li: li[i.upper()]+=1
		else: li['*']+=1
		if i == 'n' or i == 'N': Nc += 1
	
	
	avn+=1
	if args.verbose: print seq_record.id,  len(seq_record.seq), 'bps (of which '+str(Nc)+' Ns -> '+str(float(Nc)/float(len(seq_record.seq))*100)+' % )' 

	
if args.snpsFile:
	for seq_record1,seq_record2 in zip(SeqIO.parse(args.fastaFile, args.type),SeqIO.parse(args.snpsFile, args.type)):
		print seq_record1.id + ' x ' + seq_record2.id
		dist = 0
		pox  = 0
		for z,w in zip(str(seq_record1.seq).upper(),str(seq_record2.seq).upper()):
			pox+=1
			if z != w:
			    dist+=1
			    if args.verbose:print "\tSNP detected @position ",pox," (",z," vs ",w,")"
			
		print str(dist),' snps\n\n'
		
if args.trimlen:
	for seq_record in SeqIO.parse(args.fastaFile, args.type):
		if len(seq_record.seq) == int(args.trimlen):
			print '>'+str(seq_record.id)+' '+str(seq_record.description)
			print str(seq_record.seq)
if args.verbose: print "total:",total,'bps on',avn,'elements'

if args.tabular:
	#li2=[x/total for x in li.]
	print 'A\tT\tC\tG\tN\t-\t*\ttotal_bases'
	print str(round(float(li['A'])/float(total)*100.0,2))+'\t'+ str(round(float(li['T'])/float(total)*100.0,2))  +'\t'+  str(round(float(li['G'])/float(total)*100.0,2))  +'\t'+ str(round(float(li['C'])/float(total)*100.0,2))  +'\t'+str(round(float(li['N'])/float(total)*100.0,2))  +'\t'+ str(round(float(li['-'])/float(total)*100.0,2))+'\t'+  str(round(float(li['*'])/float(total)*100.0,2)) +'\t'+str(total)

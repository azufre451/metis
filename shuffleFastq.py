from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import argparse,sys,random

parser = argparse.ArgumentParser(description="FASTX files shuffler")
parser.add_argument("input", help="file containing the sequences")
parser.add_argument("output", help="file outputing the sequences")
parser.add_argument("--type", default="fasta", help="input file type (fasta|fastq)")
args=parser.parse_args()
 
seqList = [] 

i =0
for seq_record in SeqIO.parse(args.input, args.type):
	seqList.append(seq_record) 
	i+=1
	if i % 100000 == 0: print i,"seqs done"
 
print "sequences read.",len(seqList)
random.shuffle(seqList)
SeqIO.write(seqList,args.output,args.type) 

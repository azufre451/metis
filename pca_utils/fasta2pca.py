#!/usr/bin/env python

from Bio import SeqIO
import argparse,sys


parser = argparse.ArgumentParser(description='Converts a FASTA or FASTQ sequence file in matrix of binary values')

parser.add_argument("-i", required=True, help="FASTA file containing the sequences, stdin if - is specified")
parser.add_argument("-o", required=True, help="Output files, stdout if - is specified")
parser.add_argument("--type", default="fasta", help="input file type (fasta|fastq)") 
args=parser.parse_args()

with (open(args.o,'w') if args.o != '-' else sys.stdout) as outstream:
	for seq_record in SeqIO.parse(args.i if args.i != '-' else sys.stdin, args.type):
		nseq='\t'.join([a for a in str(seq_record.seq).replace('A','1000').replace('C','0100').replace('G','0010').replace('T','0001')])
		outstream.write(str(seq_record.id)+'\t'+nseq+'\n')

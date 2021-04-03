#!/usr/bin/env python

from Bio import SeqIO
#import sys
import argparse
import bz2



parser = argparse.ArgumentParser(description='Parses a FASTA file, reports seq IDs and lengths')
parser.add_argument('ifile',help="the input FASTA file (can be .bz2)")

args = parser.parse_args()

if(args.ifile.endswith('.bz2')):
	hdl=bz2.open(args.ifile,'rt')
else:
	hdl=open(args.ifile,'r')

for seq in SeqIO.parse(hdl,'fasta'):
	print(seq.id,len(seq.seq))

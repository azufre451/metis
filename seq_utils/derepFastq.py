#!/usr/bin/env python3

import sys
import os
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='Dereplicate FASTQ file from STDIN')
parser.add_argument('outfile', help='output file')
args = parser.parse_args()

seen={}
utp=[]
with open(args.outfile, "a") as output_handle:
	for read in SeqIO.parse(sys.stdin,'fastq'):
		if read.id not in seen:
			seen[read.id] = 1
			utp.append(read)
			#print(read.id,'NEW',read.seq)


		if len(utp) == 10000:
			#print(len(utp))
			SeqIO.write(utp,output_handle,'fastq')
			utp=[]

	if len(utp) > 0:
		SeqIO.write(utp,output_handle,'fastq')

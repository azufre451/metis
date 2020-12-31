#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import argparse,sys
import argparse
import os,bz2

#@HISEQ:153:C8D3VACXX:1:1304:15968:33760 1:N:0:CTTACTAGGTAAGGAG


parser = argparse.ArgumentParser(description="FASTQ Paired-End files separation based on ILLUMINA HISEQ reads signatures")

parser.add_argument("input", help="in")
parser.add_argument("prefix", help="prefix for .R1 .R2 .UN files")
parser.add_argument('-d',action='store_true')
args = parser.parse_args()

r = {}

#for record in SeqIO.parse(sys.stdin,'fastq'):
liner=[]
ptl=0

if args.d:
	fil = bz2.BZ2File(args.input,'r')
	filR1 = bz2.BZ2File(args.prefix+'.R1.fastq.bz2','w')
	filR2 = bz2.BZ2File(args.prefix+'.R2.fastq.bz2','w')
	filUP = bz2.BZ2File(args.prefix+'.UP.fastq.bz2','w')
else:
	fil = open(args.input,'r')
	filR1 = open(args.prefix+'.R1.fastq','w')
	filR2 = open(args.prefix+'.R2.fastq','w')
	filUP = open(args.prefix+'.UP.fastq','w')

for line in fil:

	if ptl == 4:
		read_name = liner[0]
		read_id = read_name.split(' ')[0]
		#print line.strip()
		flavour = liner[0].split(' ')[1].split(':')[0]  #TODO optimize for other cases, now only supports HISEQ:153:C8D3VACXX:1:1304:15968:33760 1:N:0:CTTACTAGGTAAGGAG
		if read_id not in r: r[read_id] = {}
		r[read_id][flavour] = read_name
		ptl=0
		liner=[]
	
	liner.append(line.strip()) 
	ptl += 1


R1s = set([rec['1'] for read,rec in r.items() if len(rec) == 2])
R2s = set([rec['2'] for read,rec in r.items() if len(rec) == 2])
UPs = set([rec.values()[0] for read,rec in r.items() if len(rec) == 1])

del r
print "OK1"
fil.seek(0)


#read of group 1
 

ptl=0
liner=[]
for line in fil:

	if ptl == 4:	
		
		read_name = liner[0]
		read_id = read_name.split(' ')[0]

		if read_name in R1s: filR1.write('\n'.join(liner)+'\n')
		elif read_name in R2s: filR2.write('\n'.join(liner)+'\n')
		elif read_name in UPs: filUP.write('\n'.join(liner)+'\n')

		ptl=0
		liner=[]
	
	liner.append(line.strip()) 
	ptl += 1

fil.close()
filR1.close()
filR2.close()
filUP.close()

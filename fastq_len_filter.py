#!/usr/bin/env python

from Bio import SeqIO
import argparse as ap
import sys
import numpy as np 

def read_params(args):
	p = ap.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
	p.add_argument('--min_len', required = True, default = None, type = int)
	p.add_argument('--min_qual', default = 0, type = int)
	p.add_argument('--no_anonim', action="store_true")
	p.add_argument('--count')

	return vars(p.parse_args())

screenQual=range(20,31)
qualCounter = dict((k,0) for k in screenQual)

counter = 0
allCounter = 0
rpl=[]
if __name__ == '__main__':
	args = read_params(sys.argv)
	min_len = args['min_len']
	with sys.stdout as outf:
		for r in SeqIO.parse(sys.stdin, "fastq"):
			#print r.id
			#print r.seq
			#print r.letter_annotations['phred_quality']
			#print '---'

			avQual= np.mean(r.letter_annotations['phred_quality'])
			
			if avQual >= args['min_qual'] and len(r) >= min_len:
#				print avQual	 
				if not args['no_anonim']:
					r.id = r.id+'_'+str(allCounter)

				counter+=1
				rpl.append(r)

				for qu in screenQual:
					if avQual >= qu:
						qualCounter[qu]+=1

				
				if len(rpl) % 30000 == 0:
			
					SeqIO.write(rpl, outf, "fastq")
					rpl=[]

			allCounter+=1


		if len(rpl) > 0:
			SeqIO.write(rpl, outf, "fastq")
			rpl=[]
			 
		if args['count']:
			outCount = open(args['count'],'w')
			outCount.write(str(counter)+'\t'+str(allCounter)+'\t'+'\t'.join([str(ke)+':'+str(val) for ke,val in qualCounter.items()]))
			outCount.close()


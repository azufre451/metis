import argparse
import sys

def difference(a,b):
	#hamming diff
	return sum([(1 if k!=j else 0) for k,j in zip(a,b)])

parser = argparse.ArgumentParser(description='Demultiplex FASTQ file from STDIN to STDOUT using a given barcode on the description part of FASTQ seq') 
parser.add_argument('--bcd', help='barcode',required=True)

args = parser.parse_args()


from Bio import SeqIO

tlp = []
for seq_record in SeqIO.parse(sys.stdin, "fastq"):
	rn = seq_record.id
	bcd_all = seq_record.description.split(' ')[1].split(':')[3]
	bcd_1 = bcd_all[0:8]
	bcd_2 = bcd_all[8::]

	argsbcd1 = args.bcd[0:8]
	argsbcd2 = args.bcd[8::]

	
	if difference(bcd_1,argsbcd1) <= 1 and difference(bcd_2,argsbcd2) <= 1: 
		if len(tlp) > 300000:
			SeqIO.write(tlp,sys.stdout,'fastq')
			tlp=[]
		else: 
			tlp.append(seq_record)